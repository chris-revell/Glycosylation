#%%
# flux_νₑ = (diffusive_flux_ν + advective_flux_ν)
# flux_νₑ = K₂*K₄.*Pν*∇ₑ*cᵥ - β*Pν*cₑ    where cᵥ is concentration over vertices, cₑ is concentration over edges 
# cₑ = Aᵤₚ*cᵥ
# flux_νₑ = (K₂*K₄.*Pν*∇ₑ - β*Pν*Aᵤₚ)*cᵥ
# flux_xyₑ = Dₑ*hₑ*diffusive_flux_xy
# flux_xyₑ = Dₑ*hₑ*Pxy*∇ₑ*cᵥ
# ċ = aE∇⋅flux_νₑ + a∇⋅flux_xyₑ
# ċ = a*E*∇⋅(K₂*K₄.*Pν*∇ₑ*cᵥ - β*Pν*Aᵤₚ*cᵥ) + a∇⋅(Dₑ*hₑ*Pxy*∇ₑ*cᵥ)
# Dₑ constant over edges 
# ċ = a*(E*∇⋅(K₂*K₄.*Pν*∇ₑ - β*Pν*Aᵤₚ) + 𝓓.*∇⋅(hₑ*Pxy*∇ₑ))*cᵥ

# L = -W⁻¹*Aᵀ*𝓓*l⁻¹*A .+ W⁻¹*Aᵀ*V*Aᵤₚ # Express model as a matrix operator 


# Cνν = W⁻¹*Aᵀ*Pν*l⁻¹*A
# Cν = Aᵀ*l⁻¹*Pν*Aᵤₚ
# flux_νₑ = (diffusive_flux_ν + advective_flux_ν)
# flux_νₑ = K₂*K₄.*Pν*∇ₑ*cᵥ - β*Pν*cₑ    where cᵥ is concentration over vertices, cₑ is concentration over edges 
# cₑ = Aᵤₚ*cᵥ
# flux_νₑ = (K₂*K₄.*Pν*∇ₑ - β*Pν*Aᵤₚ)*cᵥ
# flux_xyₑ = Dₑ*hₑ*diffusive_flux_xy
# flux_xyₑ = Dₑ*hₑ*Pxy*∇ₑ*cᵥ
# ċ = aE∇⋅flux_νₑ + a∇⋅flux_xyₑ
# ċ = a*E*∇⋅(K₂*K₄.*Pν*∇ₑ*cᵥ - β*Pν*Aᵤₚ*cᵥ) + a∇⋅(Dₑ*hₑ*Pxy*∇ₑ*cᵥ)
# Dₑ constant over edges 
# ċ = a*(E*∇⋅(K₂*K₄.*Pν*∇ₑ - β*Pν*Aᵤₚ) + 𝓓.*∇⋅(hₑ*Pxy*∇ₑ))*cᵥ


#


using OrdinaryDiffEq
using SparseArrays
using UnPack
using CairoMakie 
using FromFile
using DrWatson
using Printf
using SciMLOperators
using DataFrames
using Dates
using Statistics

@from "$(srcdir("MakeIncidenceMatrix.jl"))" using MakeIncidenceMatrix
@from "$(srcdir("MakeWeightMatrices.jl"))" using MakeWeightMatrices
@from "$(srcdir("Visualise.jl"))" using Visualise
@from "$(srcdir("UsefulFunctions.jl"))" using UsefulFunctions
@from "$(srcdir("CisternaWidth.jl"))" using CisternaWidth

cisternaSeriesID = 1
nSpatialDims = 2

Ωperp = 100.0  # Lumen footprint area
N     = 100         # Maximum polymer length 
k_Cd = 0.9 # Complex desorption rate
k_Ca = 1.1 # Complex adsorption rate
k_Sd = 0.9 # Substrate desorption rate
k_Sa = 1.1 # Substrate adsorption rate
k₁   = 1.0   # Complex formation forward reaction rate 
k₂   = 0.6   # Complex dissociation reverse reaction rate 
k₃   = 1.1   # Product formation
k₄   = 0.6  # Product dissociation 
E_0 = 1.0
𝓒 = 100.0
𝓢 = 100.0
D_C  = 1.0  # Monomer/polymer diffusivity
D_S  = 1.0  # Substrate diffusivity
Tᵣstar  = 100.0  # Release time
ϕ = 0.5

# PDE discretisation parameters 
Nν       = 101             # Number of discretisation points in polymerisation space
Nx       = 101             # Number of discretisation points in space
Ny       = 101             # Number of discretisation points in space
Nghost   = 1           # Number of ghost points on each side of the domain 
Nνplus   = Nν+2*Nghost # Number of discretised points including ghost points 
Nxplus   = Nx+2*Nghost # Number of discretised points including ghost points 
Nyplus   = Ny+2*Nghost # Number of discretised points including ghost points 
dimsReal = [Nν, Nx, Ny]
dimsPlus = [Nνplus, Nxplus, Nyplus]

# Generate xMax and width profile from data or function 
# xMax, mat_h = hFromData(dimsPlus; cisternaSeriesID=1)
xMax, mat_h, xs = hFromFunction(dimsPlus)
dx   = xs[2]-xs[1]
yMax = xMax
ys   = collect(range(0.0, yMax, Nyplus))
dy   = ys[2]-ys[1]
νMax = 1.0
νs   = collect(range(0.0, νMax, Nνplus)) # Positions of discretised vertices in polymerisation space 
dν   = νs[2]-νs[1]
spacing  = [dν, dx, dy]

h₀ = mean(selectdim(mat_h, 2, 1:Nxplus))

𝓔    = 2*Ωperp*E_0   # Total enzyme mass
K₃  = k₃/k₁    # Non-dimensionalised product formation rate
K₄  = k₄/k₁    # Non-dimensionalised prodict dissociation rate
δ_C = π*D_C/(k₁*𝓔)
δ_S = π*D_S/(k₁*𝓔)
Tᵣ  = k₁*𝓔*Tᵣstar/(2*Ωperp)
Ω     = h₀*Ωperp         # Lumen volume
α_C = (k_Cd*Ω)/(2*k_Ca*Ωperp) # Balance of complex in bulk to complex on membrane       units of m²?
α_S = (k_Sd*Ω)/(2*k_Sa*Ωperp) # Balance of substrate in bulk to substrate on membrane   units of m²?
C_b  = 𝓒/Ω 
S_b  = 𝓢/Ω 
C_0 = C_b*h₀/(2*(1+α_C))      # Early surface monomer concentration
S_0 = S_b*h₀/(2*(1+α_S))      # Early surface substrate concentration 
K₂  = k₂/(k₁*C_0)              # (k₂/(k₁*C_b))*((2*k_Ca*Ωperp + k_Cd*Ω)/(k_Ca*Ω)) # Non-dimensionalised complex formation net reaction rate
σ   = (k_Sa*S_b*(2*k_Ca*Ωperp + k_Cd*Ω)) / (k_Ca*C_b*(2*k_Sa*Ωperp + k_Sd*Ω))
ϵ   = 𝓔*(2*k_Ca*Ωperp + k_Cd*Ω) / (2*k_Ca*C_b*Ωperp)
𝓓   = α_C*δ_C*N^2*(K₂ + σ*K₃)
β = N*(σ*K₃ - K₂*K₄)

λ = (𝓢/(2*Ωperp))*(k₁*k₃/(k₂*k₄))
@show λ

#%%

# Create directory for run data labelled with current time.
paramsName = @savename cisternaSeriesID K₂ K₃ K₄ α_C δ_C σ N Tᵣ nSpatialDims
folderName = "$(paramsName)_$(Dates.format(Dates.now(),"yy-mm-dd-HH-MM-SS"))"
# Create frames subdirectory to store system state at each output time
subFolder = ""
mkpath(datadir("sims",subFolder,folderName))

#%%
A   = makeIncidenceMatrix3D(dimsPlus)
Ā   = abs.(A)
Aᵀ  = transpose(A)
Aᵤₚ = dropzeros((Ā-A).÷2)

nVerts  = Nνplus*Nxplus*Nyplus      # Total number of vertices 
nEdgesi = (Nνplus-1)*Nxplus*Nyplus  # Number of i-directed edges (ν, in this case)
nEdgesj = Nνplus*(Nxplus-1)*Nyplus  # Number of j-directed edges (x, in this case)
nEdgesk = Nνplus*Nxplus*(Nyplus-1)  # Number of j-directed edges (x, in this case)
nEdges  = nEdgesi+nEdgesj+nEdgesk   # Total number of edges over all dimensions 

# Ghost point masks
ghostVertexMaskVec = makeGhostVertexMask(dimsPlus)
ghostVertexMaskSparse = spdiagm(ghostVertexMaskVec)
ghostEdgeMaskVec = makeGhostEdgeMask(dimsPlus)
ghostEdgeMaskSparse = spdiagm(ghostEdgeMaskVec)

# Matrices for picking out ν and xy directions in derivatives 
Pν  = ghostEdgeMaskSparse*spdiagm(vcat(ones(Int64, nEdgesi), zeros(Int64, nEdgesj), zeros(Int64, nEdgesk)))   # Diagonal sparse matrix to exclude all xy edges 
Pxy  = ghostEdgeMaskSparse*spdiagm(vcat(zeros(Int64, nEdgesi), ones(Int64, nEdgesj), ones(Int64, nEdgesk)))   # Diagonal sparse matrix to exclude all ν edges 

# Weights
W = vertexVolumeWeightsMatrix(dimsPlus, spacing)
W⁻¹ =  vertexVolumeWeightsInverseMatrix(dimsPlus, spacing)
l⁻¹ = edgeLengthInverseMatrix(dimsPlus, spacing)

∇ₑ = l⁻¹*A       # Gradient operator giving gradient on each edge
∇cdot = -W⁻¹*Aᵀ  # Divergence operator giving divergence on each vertex calculated from edges 

# Diffusivity field over edges 
Aperpₑ = edgePerpendicularAreaMatrix(dimsPlus, spacing) # Diagonal matrix of areas perpendicular to each edge, meaning the area through which diffusive flux in the direction of a given edge passes
# Set no-flux boundary conditions by enforcing zero diffusivity in edges connection ghost points
𝓓ₑ       = 𝓓.*ghostEdgeMaskSparse*Aperpₑ # Sparse diagonal matrix of diffusivities over edges 

# Diagonal matrices of compartment thickness h over all vertices hᵥ
# Also diagonal matrix of thickness over edges, formed by taking mean of h at adjacent vertices 0.5.*Ā*hᵥ
hᵥ_vec = reshape(mat_h, nVerts)         # Cisternal thickness evaluated over vertices 
hₑ_vec = 0.5.*Ā*hᵥ_vec                  # Cisternal thickness evaluated over edges (mean of adjacent vertices)
hᵥ = spdiagm(hᵥ_vec)                    # Cisternal thickness over vertices, as a sparse diagonal matrix
hₑ = spdiagm(hₑ_vec)                    # Cisternal thickness over edges, as a sparse diagonal matrix
aᵥ = spdiagm(1.0./(1.0 .+ α_C.*hᵥ_vec)) # Prefactor 1/(1+α_C*hᵥ(x)) evaluated over vertices, packaged into a sparse diagonal matrix for convenience
aₑ = spdiagm(1.0./(1.0 .+ α_C.*hₑ_vec)) # Prefactor 1/(1+α_C*hₑ(x)) evaluated over edges, packaged into a sparse diagonal matrix for convenience

# Initial conditions using Gaussian
u0fun(x, μx, σx, y, μy, σy, z, μz, σz) = exp(-(x-μx)^2/σx^2 - (y-μy)^2/σy^2 - (z-μz)^2/σz^2)
μνu0 = 0.0; σνu0 = νMax/20.0
μxu0 = xMax/2.0; σxu0 = 10.0*xMax
μyu0 = xMax/2.0; σyu0 = 10.0*xMax
uMat = zeros(Float64, dimsPlus...)
for yy=1:Nyplus, xx=1:Nxplus, νν=1:Nνplus
    uMat[νν, xx, yy] = u0fun(νs[νν], μνu0, σνu0, xs[xx], μxu0, σxu0, ys[yy], μyu0, σyu0)
end
u0 = reshape(uMat, nVerts)
u0[ghostVertexMaskVec.!=true] .= 0.0
# For integration to normalise the number of monomers, we need to multiply the concentration at each point by the ν value of that point
νMat = ones(dimsPlus...)
for ii=1:dimsPlus[1]
    selectdim(νMat, 1, ii) .*= (ii-1)
    # νMat[ii,:].*=(ii-1)
end
νSparse = spdiagm(reshape(νMat, nVerts))
integ = sum(ghostVertexMaskSparse*W*νSparse*u0)
u0 .*= 𝓒/integ

# Set value of Fₑ at each point in space
matFₑ = ones(Float64, dimsPlus[2:end]...)
integF = prod(spacing[2:end])*sum(matFₑ[2:end-1,2:end-1])
# Ensure integral of Fₑ over space is π
matFₑ .*= π/integF
matE = zeros(dimsPlus[2:end])
Esparse = spzeros(nVerts, nVerts)
E!(u0, dimsPlus, Esparse, matE, matFₑ, K₂, dν)

# PDE operator components
L1 = aᵥ*∇cdot*(K₂*K₄.*Pν*∇ₑ - β.*Pν*Aᵤₚ)
L2 = aᵥ*∇cdot*(hₑ*Pxy*𝓓ₑ*∇ₑ)
p = (L1=L1, L2=L2, u0=u0, dimsPlus=dimsPlus, Esparse=Esparse, matE=matE, matFₑ=matFₑ, K₂=K₂, dν=dν)
L = MatrixOperator(Esparse*L1.+L2, update_func! = updateOperator!)
prob = ODEProblem(L, u0, (0.0, Tᵣ), p)
println("starting solver")
sol = solve(prob, Vern9(), saveat=Tᵣ/100.0)

#%%

println("finished sim")

# productionHeatmap3D(ϕ, sol.u, sol.t, xs, νs, dimsReal, ghostVertexMaskVec, W; subFolder=subFolder, folderName=folderName)
spaceIntegralOver_ν_Movie(sol.u, sol.t, xs, νs, dimsReal, Nghost, W, ghostVertexMaskVec; subFolder=subFolder, folderName=folderName)