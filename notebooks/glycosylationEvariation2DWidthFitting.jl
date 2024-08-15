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

@from "$(srcdir("MakeIncidenceMatrix.jl"))" using MakeIncidenceMatrix
@from "$(srcdir("MakeWeightMatrices.jl"))" using MakeWeightMatrices
@from "$(srcdir("Visualise.jl"))" using Visualise
@from "$(srcdir("UsefulFunctions.jl"))" using UsefulFunctions
@from "$(srcdir("CisternaWidth.jl"))" using CisternaWidth

cisternaSeriesID = 1

# PDE discretisation parameters 
Nx       = 101             # Number of discretisation points in space
Nν       = 101             # Number of discretisation points in polymerisation space
Nghost   = 1           # Number of ghost points on each side of the domain 
Nνplus   = Nν+2*Nghost # Number of discretised points including ghost points 
Nxplus   = Nx+2*Nghost # Number of discretised points including ghost points 
dimsPlus = (Nνplus, Nxplus)

# Generate xMax and width profile from data or function 
# xMax, mat_h = hFromData(dimsPlus; cisternaSeriesID=1)
xMax, mat_h, xs = hFromFunction(dimsPlus)
dx   = xs[2]-xs[1]
νMax = 1.0
νs   = collect(range(0.0, νMax, Nνplus)) # Positions of discretised vertices in polymerisation space 
dν   = νs[2]-νs[1]
spacing  = (dν, dx)

K₁ = 1.0
K₂ = 1.0
K₃ = 2.0
K₄ = 1.0  
α_C = 100.0
δ_C = 1.0
σ = 10.0
N = 100
β = N*(σ*K₃ - K₂*K₄)
𝓓 = α_C*δ_C*N^2*(K₂+σ*K₃)
Tᵣ = 0.2

#%%

# Create directory for run data labelled with current time.
paramsName = @savename cisternaSeriesID K₂ K₃ K₄ α_C δ_C σ N Tᵣ
folderName = "$(paramsName)_$(Dates.format(Dates.now(),"yy-mm-dd-HH-MM-SS"))"
# Create frames subdirectory to store system state at each output time
subFolder = ""
mkpath(datadir("sims",subFolder,folderName))

u0fun(x, μx, σx, y, μy, σy) = exp(-(x-μx)^2/σx^2 - (y-μy)^2/σy^2)
μνu0 = 0.0; σνu0 = νMax/10.0
μxu0 = xMax/2.0; σxu0 = 10.0*xMax

fFun(x, μx, σx) = 0.1 #+ exp(-(x-μx)^2/σx^2)
μxF = xMax/2.0; σxF=xMax/10.0

#%%

A   = makeIncidenceMatrix3D(Nνplus, Nxplus, 1)
Ā   = abs.(A)
Aᵀ  = transpose(A)
Aᵤₚ = dropzeros((Ā-A).÷2)

nVerts  = Nνplus*Nxplus       # Total number of vertices 
nEdgesi = (Nνplus-1)*Nxplus  # Number of i-directed edges (ν, in this case)
nEdgesj = Nνplus*(Nxplus-1)  # Number of j-directed edges (x, in this case)
nEdges  = nEdgesi+nEdgesj     # Total number of edges over all dimensions 

# Ghost point masks
ghostVertexMask = makeGhostVertexMask(dimsPlus)
ghostVertexMaskSparse = spdiagm(ghostVertexMask)
ghostEdgeMask = makeGhostEdgeMask(dimsPlus)
ghostEdgeMaskSparse = spdiagm(ghostEdgeMask)

# Weights
W = vertexVolumeWeightsMatrix(dimsPlus, spacing)
W⁻¹ =  vertexVolumeWeightsInverseMatrix(dimsPlus, spacing)
l⁻¹ = edgeLengthInverseMatrix(dimsPlus, spacing)

# Diagonal matrices of compartment thickness h over all vertices hᵥ
# Also diagonal matrix of thickness over edges, formed by taking mean of h at adjacent vertices 0.5.*Ā*hᵥ
hᵥ_vec = reshape(mat_h, nVerts)         # Cisternal thickness evaluated over vertices 
hₑ_vec = 0.5.*Ā*hᵥ_vec                  # Cisternal thickness evaluated over edges (mean of adjacent vertices)
hᵥ = spdiagm(hᵥ_vec)                    # Cisternal thickness over vertices, as a sparse diagonal matrix
hₑ = spdiagm(hₑ_vec)                    # Cisternal thickness over edges, as a sparse diagonal matrix
aᵥ = spdiagm(1.0./(1.0 .+ α_C.*hᵥ_vec)) # Prefactor 1/(1+α_C*hᵥ(x)) evaluated over vertices, packaged into a sparse diagonal matrix for convenience
aₑ = spdiagm(1.0./(1.0 .+ α_C.*hₑ_vec)) # Prefactor 1/(1+α_C*hₑ(x)) evaluated over edges, packaged into a sparse diagonal matrix for convenience

# Velocity field 
V_i = fill(β, (Nνplus-1, Nxplus))
V_j = fill(0.0, (Nνplus, Nxplus-1))
Vvec = vcat(reshape(V_i, nEdgesi), reshape(V_j, nEdgesj))
V = ghostEdgeMaskSparse*spdiagm(Vvec)*aₑ   # Diagonal matrix of advection velocities at each edge

# Diffusivity field over edges 
# Set no-flux boundary conditions by enforcing zero diffusivity in edges connection ghost points
𝓓_i  = fill(dx*K₂*K₄, (Nνplus-1, Nxplus))
𝓓_j  = fill(dν*K₂*K₄, (Nνplus, Nxplus-1))
𝓓vec = vcat(reshape(𝓓_i, nEdgesi), reshape(𝓓_j, nEdgesj))
𝓓    = ghostEdgeMaskSparse*spdiagm(𝓓vec)*aₑ # Diagonal matrix of advection velocities at each edge

# Matrices for picking out ν and xy directions in derivatives 
P  = ghostEdgeMaskSparse*spdiagm(vcat(ones(Int64, nEdgesi), ones(Int64, nEdgesj)))     # Diagonal sparse matrix to exclude all edges adjacent to ghost points  
Pν = ghostEdgeMaskSparse*spdiagm(vcat(ones(Int64, nEdgesi), zeros(Int64, nEdgesj)))   # Diagonal sparse matrix to exclude all xy edges and ν edges adjacent to ghost points  
Px = ghostEdgeMaskSparse*spdiagm(vcat(zeros(Int64, nEdgesi), ones(Int64, nEdgesj)))   # Diagonal sparse matrix to exclude all ν edges and xy edges adjacent to ghost points 

# Diagonal matrix of edge lengths
l_i = fill(dν, (Nνplus-1, Nxplus))
l_j = fill(dx, (Nνplus, Nxplus-1))
lvec = vcat(reshape(l_i, nEdgesi), reshape(l_j, nEdgesj))
l = spdiagm(lvec)
l⁻¹ = spdiagm(1.0./lvec)

# Initial conditions using Gaussian
uMat = zeros(Float64, Nνplus, Nxplus)
for xx=1:Nxplus, νν=1:Nνplus
    uMat[νν, xx] = u0fun(νs[νν], μνu0, σνu0, xs[xx], μxu0, σxu0)            
end
u0 = reshape(uMat, nVerts)
u0[ghostVertexMask.!=true] .= 0.0
integ = sum(W*u0)
u0 ./= integ

∇ₑ = l⁻¹*A       # Gradient operator giving gradient on each edge
∇cdot = -W⁻¹*Aᵀ  # Divergence operator giving divergence on each vertex calculated from edges 

matFₑ = zeros(Nνplus, Nxplus)
for j=1:Nxplus
    matFₑ[:, j] .= fFun(xs[j], μxF, σxF)
    # selectdim(matFₑ, 2, j) .= fFun(xs[j], μxF, σxF)
    # matFₑ[i] = 1.0
end
matE = zeros(Nνplus, Nxplus)
E = spdiagm(reshape(matE, nVerts))

L1 = aᵥ*∇cdot*(K₂*K₄.*Pν*∇ₑ - β.*Pν*Aᵤₚ)
L2 = aᵥ*∇cdot*(𝓓*hₑ*Px*∇ₑ)

p = (L1 = L1,
    L2 = L2,
    Nνplus = Nνplus,
    Nxplus = Nxplus,
    K₂ = K₂,
    matE = matE,
    E = E,
    matFₑ = matFₑ)

function update_func!(L, u, p, t)
    @unpack L1,
        L2,
        Nνplus,
        Nxplus,
        K₂,
        matE,
        E,
        matFₑ = p

    cs = reshape(u, (Nνplus, Nxplus))     
    for j = 1:Nxplus
        integrationFactor = K₂/(K₂ + simpsonsRule(cs[:,j]))
        matE[:,j] .= matFₑ[:,j].*integrationFactor
    end
    E .= spdiagm(reshape(matE, nVerts)) 
    L .= E*L1 .+ L2
end

L = MatrixOperator(E*L1.+L2, update_func! = update_func!)
prob = ODEProblem(L, u0, (0.0, Tᵣ), p)
sol = solve(prob, Vern9(), saveat=Tᵣ/100.0)

#%%

concentrationSurfaceMovie(sol.u, sol.t, xs, νs, (Nν,Nx), Nghost, ghostVertexMask; subFolder=subFolder, folderName=folderName)

spaceIntegralOver_ν_Movie(sol.u, sol.t, xs, νs, (Nν,Nx), Nghost, W, ghostVertexMask; subFolder=subFolder, folderName=folderName)



# N = 100.0

# σ = 1000.0 # >> Non-dimensionalised substrate availability 
# α_C = 100.0 # >> 1 measures the solubility of the cargo; equivalently, inverse of affinity for the membrane. α_C >> 1 suggests cargo adsorbs weakly onto membrane

# K₃ = 1.1/(N*σ) # Non-dimensionalised complex dissociation Q->C+E rate relative to complex formation C+E->Q rate (k₁ -> K₁=1.0) 
# K₂ = 0.9999/N # Non-dimensionalised product formation Q+S->C+E rate relative to complex formation C+E->Q rate (k₁ -> K₁=1.0)
# K₄ = 0.9999/N # Non-dimensionalised product dissociation C+E-> Q+S rate relative to complex formation C+E->Q rate (k₁ -> K₁=1.0)
# β  = N*(σ*K₃ - K₂*K₄) # ≈ 1.0 Balance of production and dissociation of cargo

# δ_C = 0.00001

# # 𝓓  = 1.0 # α_C*δ_C*N^2*(K₂ + σ*K₃) Nondimensionalised diffusion constant 
# 𝓓  = α_C*δ_C*N^2*(K₂ + σ*K₃) # Nondimensionalised diffusion constant 

# Tᵣ = 1.0

