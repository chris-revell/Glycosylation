

dimsPlus = (10, 12, 14)
spacing = (0.1,0.2,0.3)

# Incidence matrices 
A   = makeIncidenceMatrix3D(dimsPlus...)
Ā   = abs.(A)
Aᵀ  = transpose(A)
Aᵤₚ = dropzeros((Ā-A).÷2)

# Number of vertices and number of edges, total and in each dimension
nVerts  = prod(dimsPlus)      # Total number of vertices 

nEdgesijk = Int64[]
for i=1:length(dimsPlus)
    push!(nEdgesijk, (dimsPlus[i]-1)*prod(dimsPlus[1:end .!= i]))
end
nEdges  = sum(nEdgesijk)     # Total number of edges over all dimensions 

# Ghost point masks; vectors and sparse diagonal matrices to exclude ghost points and edges connected to ghost points 
ghostVertexMaskVec = makeGhostVertexMask(dimsPlus)
ghostVertexMaskSparse = spdiagm(ghostVertexMaskVec)
ghostEdgeMaskVec = makeGhostEdgeMask(dimsPlus)
ghostEdgeMaskSparse = spdiagm(ghostEdgeMaskVec)

# Weights
W   = vertexVolumeWeightsMatrix(dimsPlus, spacing)
W⁻¹ = vertexVolumeWeightsInverseMatrix(dimsPlus, spacing)
l⁻¹ = edgeLengthInverseMatrix(dimsPlus, spacing) # Diagonal matrix of edge lengths

# Gradient operators 
∇ₑ = l⁻¹*A       # Gradient operator giving gradient on each edge
∇cdot = -W⁻¹*Aᵀ  # Divergence operator giving divergence on each vertex calculated from edges 

# Diagonal matrices of compartment thickness h over all vertices hᵥ
# Also diagonal matrix of thickness over edges, formed by taking mean of h at adjacent vertices 0.5.*Ā*hᵥ
hᵥ_vec = reshape(mat_h, nVerts)         # Cisternal thickness evaluated over vertices 
hₑ_vec = 0.5.*Ā*hᵥ_vec                  # Cisternal thickness evaluated over edges (mean of adjacent vertices)
hᵥ = spdiagm(hᵥ_vec)                    # Cisternal thickness over vertices, as a sparse diagonal matrix
hₑ = spdiagm(hₑ_vec)                    # Cisternal thickness over edges, as a sparse diagonal matrix
aᵥ = spdiagm(1.0./(1.0 .+ α_C.*hᵥ_vec)) # Prefactor 1/(1+α_C*hᵥ(x)) evaluated over vertices, packaged into a sparse diagonal matrix for convenience
aₑ = spdiagm(1.0./(1.0 .+ α_C.*hₑ_vec)) # Prefactor 1/(1+α_C*hₑ(x)) evaluated over edges, packaged into a sparse diagonal matrix for convenience

# Velocity field 
Vijk = []
sijk = (β, 0, 0)
for i=1:
V_i = fill(β, (dimsPlus[1]-1, dimsPlus[2:end]...))
Vvec = reshapew(V_i, nEdgesi), reshape(V_j, nEdgesj))
Vvec = vcat(reshape(V_i, nEdgesi), reshape(V_j, nEdgesj))
V = ghostEdgeMaskSparse*spdiagm(Vvec)*aₑ   # Diagonal matrix of advection velocities at each edge

# Diffusivity field over edges 
# Set no-flux boundary conditions by enforcing zero diffusivity in edges connection ghost points
D_i  = fill(dx*K₂*K₄, (Nνplus-1, Nxplus))
D_j  = fill(dν*K₂*K₄, (Nνplus, Nxplus-1))
Dvec = vcat(reshape(D_i, nEdgesi), reshape(D_j, nEdgesj))
D    = ghostEdgeMaskSparse*spdiagm(Dvec)*aₑ # Diagonal matrix of advection velocities at each edge

# Matrices for picking out ν and xy directions in derivatives 
# P  = ghostEdgeMaskSparse*spdiagm(ones(Int64, nEdges))     # Diagonal sparse matrix to exclude all edges adjacent to ghost points  
Pν = ghostEdgeMaskSparse*spdiagm(vcat(ones(Int64, nEdgesi), zeros(Int64, nEdgesj)))   # Diagonal sparse matrix to exclude all xy edges and ν edges adjacent to ghost points  
Px = ghostEdgeMaskSparse*spdiagm(vcat(zeros(Int64, nEdgesi), ones(Int64, nEdgesj)))   # Diagonal sparse matrix to exclude all ν edges and xy edges adjacent to ghost points 
