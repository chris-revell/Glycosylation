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

nSpatialDims = 2

# PDE discretisation parameters 
Nν       = 101             # Number of discretisation points in polymerisation space
Nx       = 101             # Number of discretisation points in space
Ny       = 101             # Number of discretisation points in space
Nghost   = 1           # Number of ghost points on each side of the domain 
Nνplus   = Nν+2*Nghost # Number of discretised points including ghost points 
Nxplus   = Nx+2*Nghost # Number of discretised points including ghost points 
Nyplus   = Ny+2*Nghost # Number of discretised points including ghost points 
dimsPlus = (Nνplus, Nxplus, Nyplus)

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
spacing  = (dν, dx, dy)

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
Tᵣ = 0.01

#%%

# Create directory for run data labelled with current time.
paramsName = @savename cisternaSeriesID K₂ K₃ K₄ α_C δ_C σ N Tᵣ nSpatialDims
folderName = "$(paramsName)_$(Dates.format(Dates.now(),"yy-mm-dd-HH-MM-SS"))"
# Create frames subdirectory to store system state at each output time
subFolder = ""
mkpath(datadir("sims",subFolder,folderName))

#%%

A   = makeIncidenceMatrix3D(Nνplus, Nxplus, Nyplus)
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

# Weights
W = vertexVolumeWeightsMatrix(dimsPlus, spacing)
W⁻¹ =  vertexVolumeWeightsInverseMatrix(dimsPlus, spacing)
l⁻¹ = edgeLengthInverseMatrix(dimsPlus, spacing)

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
V_i = fill(β, (Nνplus-1, Nxplus, Nyplus))
V_j = fill(0.0, (Nνplus, Nxplus-1, Nyplus))
V_k = fill(0.0, (Nνplus, Nxplus, Nyplus-1))
Vvec = vcat(reshape(V_i, nEdgesi), reshape(V_j, nEdgesj), reshape(V_k, nEdgesk))
V = ghostEdgeMaskSparse*spdiagm(Vvec)*aₑ   # Diagonal matrix of advection velocities at each edge

# Diffusivity field over edges 
# Set no-flux boundary conditions by enforcing zero diffusivity in edges connection ghost points
𝓓_i  = fill(dx*dy*K₂*K₄, (Nνplus-1, Nxplus, Nyplus))
𝓓_j  = fill(dν*dy*K₂*K₄, (Nνplus, Nxplus-1, Nyplus))
𝓓_k  = fill(dν*dx*K₂*K₄, (Nνplus, Nxplus, Nyplus-1))
𝓓vec = vcat(reshape(𝓓_i, nEdgesi), reshape(𝓓_j, nEdgesj), reshape(𝓓_k, nEdgesk))
𝓓    = ghostEdgeMaskSparse*spdiagm(𝓓vec)*aₑ # Diagonal matrix of advection velocities at each edge

# Matrices for picking out ν and xy directions in derivatives 
Pν = ghostEdgeMaskSparse*spdiagm(vcat(ones(Int64, nEdgesi), zeros(Int64, nEdgesj), zeros(Int64, nEdgesk)))   # Diagonal sparse matrix to exclude all xy edges and ν edges adjacent to ghost points  
Pxy = ghostEdgeMaskSparse*spdiagm(vcat(zeros(Int64, nEdgesi), ones(Int64, nEdgesj), ones(Int64, nEdgesk)))   # Diagonal sparse matrix to exclude all ν edges and xy edges adjacent to ghost points 

# Initial conditions using Gaussian
u0fun(x, μx, σx, y, μy, σy, z, μz, σz) = exp(-(x-μx)^2/σx^2 - (y-μy)^2/σy^2 - (z-μz)^2/σz^2)
μνu0 = 0.0; σνu0 = νMax/10.0
μxu0 = xMax/2.0; σxu0 = 10.0*xMax
μyu0 = xMax/2.0; σyu0 = 10.0*xMax
uMat = zeros(Float64, Nνplus, Nxplus, Nyplus)
for yy=1:Nyplus, xx=1:Nxplus, νν=1:Nνplus
    uMat[νν, xx, yy] = u0fun(νs[νν], μνu0, σνu0, xs[xx], μxu0, σxu0, ys[yy], μyu0, σyu0)
end
u0 = reshape(uMat, nVerts)
u0[ghostVertexMaskVec.!=true] .= 0.0
integ = sum(W*u0)
u0 ./= integ

matFₑ = ones(Float64, Nνplus, Nxplus, Nyplus)
vecFₑ = ones(Float64, nVerts)
integ = sum((W*vecFₑ)[ghostVertexMaskVec])
vecFₑ .*= π/integ
matE = zeros(Nνplus, Nxplus, Nyplus)
E = spdiagm(zeros(nVerts))

L1 = aᵥ*∇cdot*(K₂*K₄.*Pν*∇ₑ - β.*Pν*Aᵤₚ)
L2 = aᵥ*∇cdot*(𝓓*hₑ*Pxy*∇ₑ)

p = (L1 = L1,
    L2 = L2,
    K₂ = K₂,
    matE = matE,
    E = E,
    matFₑ = matFₑ,
    dν = dν)

function update_func!(L, u, p, t)
    @unpack L1,
        L2,
        K₂,
        matE,
        E,
        matFₑ,
        dν = p

    cs = reshape(u, size(matE))
    
    matE[1,:,:] .= matFₑ[1,:,:].*(K₂./(K₂ .+ dν.*sum(cs[2:end-1,:,:], dims=1)[1,:,:]))
    for i=2:size(matE)[1]
        matE[i, :, :] .= matE[1, :, :]
    end
    
    E .= spdiagm(reshape(matE, nVerts)) 
    L .= E*L1 .+ L2
end




p = (L1 = L1,
    L2 = L2,
    Nνplus = Nνplus,
    Nxplus = Nxplus,
    Nyplus = Nyplus,
    K₂ = K₂,
    matE = matE,
    E = E,
    matFₑ,
    νs = matFₑ)

function update_func!(L, u, p, t)
    @unpack L1,
        L2,
        Nνplus,
        Nxplus,
        Nyplus,
        K₂,
        matE,
        E,
        matFₑ,
        νs = p

    cs = reshape(u, (Nνplus, Nxplus, Nyplus))     
    for kk=1:Nyplus, jj=1:Nxplus
        integrationFactor = K₂/(K₂ + trapeziumRule(cs[:,jj,kk], νs))
        matE[:,jj,kk] .= matFₑ[:,jj,kk].*integrationFactor
    end
    E .= spdiagm(reshape(matE, nVerts))
    L .= E*L1 .+ L2
end

L = MatrixOperator(E*L1.+L2, update_func! = update_func!)
prob = ODEProblem(L, u0, (0.0, Tᵣ), p)
sol = solve(prob, Vern9(), saveat=Tᵣ/100.0)

#%%

# concentrationSurfaceMovie(sol.u, sol.t, xs, νs, (Nν,Nx), Nghost, ghostVertexMaskVec; subFolder=subFolder, folderName=folderName)

# spaceIntegralOver_ν_Movie(sol.u, sol.t, xs, νs, (Nν,Nx), Nghost, W, ghostVertexMaskVec; subFolder=subFolder, folderName=folderName)

println("finished sim")

ϕ = 0.75

productionHeatmap3D(ϕ, sol.u, sol.t, xs, νs, (Nν, Nx, Ny), ghostVertexMaskVec, W; subFolder=subFolder, folderName=folderName)


# fFun(x, μx, σx) = 1.0 #+ exp(-(x-μx)^2/σx^2)
# μxF = xMax/2.0; σxF=xMax/10.0