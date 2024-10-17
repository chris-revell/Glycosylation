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
using Dates
using Statistics
using InvertedIndices
using LinearAlgebra
@from "$(srcdir("Glycosylation.jl"))" using Glycosylation
@from "$(srcdir("Visualise.jl"))" using Visualise
@from "$(srcdir("UsefulFunctions.jl"))" using UsefulFunctions
@from "$(srcdir("MakeWeightMatrices.jl"))" using MakeWeightMatrices
@from "$(srcdir("MakeIncidenceMatrix.jl"))" using MakeIncidenceMatrix
@from "$(srcdir("DerivedParameters.jl"))" using DerivedParameters


nSpatialDims = 1

h₀ = 0.1

Ωperp = 100.0  # Lumen footprint area
N     = 100         # Maximum polymer length 
k_Cd  = 200.0 # Complex desorption rate
k_Ca  = 2.0 # Complex adsorption rate
k_Sd  = 1.0 # Substrate desorption rate
k_Sa  = 1.1 # Substrate adsorption rate
k₁    = 1.0   # Complex formation forward reaction rate 
k₂    = 0.1   # Complex dissociation reverse reaction rate 
k₃    = 1.0   # Product formation
k₄    = 1.0  # Product dissociation 
E_0   = 0.001
𝓒     = 100.0
𝓢     = 1000.0
D_C   = 0.01  # Monomer/polymer diffusivity
D_S   = 0.01  # Substrate diffusivity
Tᵣstar= 5000.0  # Release time
ϕ     = 0.5

Ngrid = 101

xMax = 100.0
xs   = collect(range(0.0, xMax, Ngrid)) # Positions of discretised vertices in space
mat_h = h₀.*ones(fill(Ngrid, nSpatialDims+1)...)

Nν   = Ngrid
Nx   = Ngrid
Ny   = Ngrid
nSpatialDims == 1 ? dims = [Nν, Nx] : dims = [Nν, Nx, Ny]
dx   = xs[2]-xs[1]
if nSpatialDims > 1 
    yMax = xMax
    ys   = collect(range(0.0, yMax, Ny))
    dy   = ys[2]-ys[1]
end
νMax = 1.0
νs   = collect(range(0.0, νMax, Nν)) # Positions of discretised vertices in polymerisation space 
dν   = νs[2]-νs[1]
nSpatialDims == 1 ? spacing  = [dν, dx] : spacing  = [dν, dx, dy]

derivedParameters(h₀, Ωperp, N, k_Cd, k_Ca, k_Sd, k_Sa, k₁, k₂, k₃, k₄, E_0, 𝓒, 𝓢, D_C, D_S, Tᵣstar; checks=true)

#%%

Nν   = Ngrid
Nx   = Ngrid
Ny   = Ngrid

nSpatialDims == 1 ? dims = [Nν, Nx] : dims = [Nν, Nx, Ny]
nSpatialDims == 1 ? dims = [Ngrid, Ngrid] : dims = [Ngrid, Ngrid, Ngrid]

# Generate xMax and width profile from data or function 
# xMax, mat_h = hFromData(dims; cisternaSeriesID=1)
xMax = maximum(xs)
dx   = xs[2]-xs[1]
if nSpatialDims > 1 
    yMax = xMax
    ys   = collect(range(0.0, yMax, Ny))
    dy   = ys[2]-ys[1]
end
νMax = 1.0
νs   = collect(range(0.0, νMax, Nν)) # Positions of discretised vertices in polymerisation space 
dν   = νs[2]-νs[1]

nSpatialDims == 1 ? spacing  = [dν, dx] : spacing  = [dν, dx, dy]

h₀ = mean(selectdim(mat_h, 1, 1))

derivedParams = derivedParameters(h₀, Ωperp, N, k_Cd, k_Ca, k_Sd, k_Sa, k₁, k₂, k₃, k₄, E_0, 𝓒, 𝓢, D_C, D_S, Tᵣstar; checks=false)
@unpack 𝓔, K₃, K₄, δ_C, δ_S, Tᵣ, Ω, α_C, α_S, C_b, S_b, C_0, S_0, K₂, σ, ϵ, 𝓓, β, K₂, L₀ = derivedParams 

λ = (𝓢/(2*Ωperp))*(k₁*k₃/(k₂*k₄))
@show λ

A   = makeIncidenceMatrix3D(dims)
Ā   = abs.(A)
Aᵀ  = transpose(A)
Aᵤₚ = dropzeros((Ā-A).÷2)

nVerts  = prod(dims)      # Total number of vertices 
dimEdgeCount = Int64[]
for i=1:length(dims)
    push!(dimEdgeCount, (dims[i]-1)*prod(dims[Not(i)]))
end
nEdges  = sum(dimEdgeCount)   # Total number of edges over all dimensions 

# Matrices for picking out ν and xy directions in derivatives 
Pν  = spdiagm(vcat(ones(Int64, dimEdgeCount[1]), zeros(Int64, sum(dimEdgeCount[2:end]))))   # Diagonal sparse matrix to exclude all xy edges 
Pxy  = spdiagm(vcat(zeros(Int64, dimEdgeCount[1]), ones(Int64, sum(dimEdgeCount[2:end]))))   # Diagonal sparse matrix to exclude all ν edges 

# Weights
W   = vertexVolumeWeightsMatrix(dims, spacing)
W⁻¹ =  vertexVolumeWeightsInverseMatrix(dims, spacing)
l⁻¹ = edgeLengthInverseMatrix(dims, spacing)

∇ₑ = l⁻¹*A       # Gradient operator giving gradient on each edge
∇cdot = -W⁻¹*Aᵀ  # Divergence operator giving divergence on each vertex calculated from edges 

# Diffusivity field over edges 
# Set no-flux boundary conditions by enforcing zero diffusivity in edges connection ghost points
Aperpₑ = edgePerpendicularAreaMatrix(dims, spacing)
𝓓ₑ     = 𝓓.*Aperpₑ # Sparse diagonal matrix of diffusivities over edges 

# Diagonal matrices of compartment thickness h over all vertices hᵥ
# Also diagonal matrix of thickness over edges, formed by taking mean of h at adjacent vertices 0.5.*Ā*hᵥ
hᵥ_vec = reshape(mat_h, nVerts)         # Cisternal thickness evaluated over vertices 
hₑ_vec = 0.5.*Ā*hᵥ_vec                  # Cisternal thickness evaluated over edges (mean of adjacent vertices)
hᵥ = spdiagm(hᵥ_vec)                    # Cisternal thickness over vertices, as a sparse diagonal matrix
hₑ = spdiagm(hₑ_vec)                    # Cisternal thickness over edges, as a sparse diagonal matrix
aᵥ = spdiagm(1.0./(1.0 .+ α_C.*hᵥ_vec)) # Prefactor 1/(1+α_C*hᵥ(x)) evaluated over vertices, packaged into a sparse diagonal matrix for convenience
# aₑ = spdiagm(1.0./(1.0 .+ α_C.*hₑ_vec)) # Prefactor 1/(1+α_C*hₑ(x)) evaluated over edges, packaged into a sparse diagonal matrix for convenience

spatialWeightsMat = reshape(diag(W))

u0fun(xs, μs, σs) = exp(-sum((xs.-μs).^2.0./σs.^2.0)) # Multidimensional Gaussian
uMat = zeros(Float64, dims...)
for ind in CartesianIndices(uMat)
    uMat[ind] = u0fun([νs[ind[1]]], [0.0], [νMax/10.0])
end
u0 = reshape(uMat, nVerts)
integ = sum(W*u0)
u0 .*= 𝓒/integ

# Set value of Fₑ at each point in space
matFₑTmp = ones(Float64, dims...)
for i=2:length(dims)
    selectdim(matFₑTmp, i, 1) .*= 0.5
    selectdim(matFₑTmp, i, dims[i]) .*= 0.5
end
integF = prod(spacing[Not(1)])*sum(selectdim(matFₑTmp, 1, 1))
# Ensure integral of Fₑ over space is π
matFₑ = (π/integF).*ones(Float64, dims[Not(1)]...)
matE = zeros(dims...)
Esparse = spzeros(nVerts, nVerts)
E!(u0, dims, Esparse, matE, matFₑ, K₂, dν)

# PDE operator components
L1 = aᵥ*∇cdot*Aperpₑ*(K₂*K₄.*Pν*∇ₑ - β.*Pν*Aᵤₚ)
L2 = aᵥ*∇cdot*(hₑ*Pxy*𝓓ₑ*∇ₑ)
p = (L1=L1, L2=L2, u0=u0, dims=dims, Esparse=Esparse, matE=matE, matFₑ=matFₑ, K₂=K₂, dν=dν)
L = MatrixOperator(Esparse*L1.+L2, update_func! = updateOperator!)

prob = ODEProblem(L, u0, (0.0, Tᵣ), p)
sol = solve(prob, Vern9(), saveat=Tᵣ/100.0, progress=true)
println("finished sim")

#%%

# Create directory for run data labelled with current time.
paramsName = @savename nSpatialDims Ωperp k_Cd k_Ca k_Sd k_Sa k₁ k₂ k₃ k₄ E_0 𝓒 𝓢 D_C D_S Tᵣstar ϕ
folderName = "$(Dates.format(Dates.now(),"yy-mm-dd-HH-MM-SS"))_$(paramsName)"
# Create frames subdirectory to store system state at each output time
subFolder = ""
mkpath(datadir("sims",subFolder,folderName))

W = vertexVolumeWeightsMatrix(dims, spacing)

if nSpatialDims==1
    concentrationSurfaceMovie(sol.u, sol.t, dims; subFolder="", folderName=folderName)
else
    spaceIntegralOver_ν_Movie(sol.u, sol.t, xs, νs, dims, W; subFolder=subFolder, folderName=folderName)
end

