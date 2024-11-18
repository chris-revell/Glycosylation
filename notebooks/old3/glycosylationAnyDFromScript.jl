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
using InvertedIndices
using XLSX
using DataFrames
using Interpolations
using Statistics
using GaussianRandomFields

@from "$(srcdir("Glycosylation.jl"))" using Glycosylation
@from "$(srcdir("Visualise.jl"))" using Visualise
@from "$(srcdir("UsefulFunctions.jl"))" using UsefulFunctions
@from "$(srcdir("MakeWeightMatrices.jl"))" using MakeWeightMatrices
@from "$(srcdir("DerivedParameters.jl"))" using DerivedParameters
@from "$(srcdir("CisternaWidth.jl"))" using CisternaWidth
@from "$(srcdir("MakeIncidenceMatrix.jl"))" using MakeIncidenceMatrix

u0Fun(xs, μs, σs) = exp(-sum((xs.-μs).^2.0./σs.^2.0)) # Multidimensional Gaussian

function hFun(dims; λ=0.1, σ=1.0)
    if length(dims) == 2
        cov = CovarianceFunction(length(dims)-1, Gaussian(λ, σ=σ))
        pts = range(0, stop=1, length=dims[2])
        grf = GaussianRandomField(cov, CirculantEmbedding(), pts, minpadding=10001)
        mat_hSlice = sample(grf)[1:dims[2]]
    else
        cov = CovarianceFunction(length(dims)-1, Gaussian(λ, σ=σ))
        pts1 = range(0, stop=1, length=dims[2])
        pts2 = range(0, stop=1, length=dims[3])
        grf = GaussianRandomField(cov, CirculantEmbedding(), pts1, pts2, minpadding=10001)
        mat_hSlice = sample(grf)[1:dims[2], 1:dims[3]]
    end
    mat_hSlice .= mat_hSlice.-mean(mat_hSlice).+1.0
    mat_h = zeros(dims...)
    for i=1:dims[1]
        selectdim(mat_h, 1, i) .= mat_hSlice
    end
    return mat_h
end

nSpatialDims = 1
Ngrid = 101
thickness = "uniform"

#%%

# Ωperp = 1000.0    # Dimensional lumen footprint area
# Ω     = 1.0      # Dimensional lumen volume 
# N     = 1000     # Maximum polymer length 
# k_Cd  = 100.0    # Dimensional complex desorption rate
# k_Ca  = 0.01     # Dimensional complex adsorption rate
# k_Sd  = 1.0      # Dimensional substrate desorption rate
# k_Sa  = 1.0      # Dimensional substrate adsorption rate
# k₁    = 1.0      # Dimensional complex formation forward reaction rate 
# k₂    = 1.0     # Dimensional complex dissociation reverse reaction rate 
# k₃    = 1.0     # Dimensional product formation
# k₄    = 1.0      # Dimensional product dissociation 
# 𝓔     = 0.001           # Dimensional total enzyme mass 
# 𝓒     = 1.0
# 𝓢     = 1000.0
# D_C   = 0.001  # Monomer/polymer diffusivity
# D_S   = 0.001  # Substrate diffusivity
# Tᵣstar= 0.01  # Release time
# ϕ     = 0.5

# dims  = fill(Ngrid, nSpatialDims+1)
# derivedParams = derivedParameters(Ω, Ωperp, N, k_Cd, k_Ca, k_Sd, k_Sa, k₁, k₂, k₃, k₄, 𝓒, 𝓢, 𝓔, D_C, D_S, Tᵣstar; checks=true)
# @unpack L₀, E₀, h₀, C_b, S_b, δ_C, δ_S, α_C, α_S, C₀, S₀, Tᵣ, K₂, K₃, K₄, σ, ϵ, 𝓓, β = derivedParams

#%%

differencing = "centre"
nSpatialDims = 1
Tᵣ = 30.0
K₂ = 1.0
K₄ = 0.0001
α_C = 1.0
𝓓 = 1.0
β = 0.1
Ngrid = 201
dims  = fill(Ngrid, nSpatialDims+1)

# PDE discretisation parameters 
nSpatialDims = length(dims)-1
    
xMax = π^(1/nSpatialDims)
xs   = collect(range(0.0, xMax, dims[2]))
dx   = xs[2]-xs[1]    
if nSpatialDims > 1 
    yMax = xMax
    ys   = collect(range(0.0, yMax, dims[3]))
    dy   = ys[2]-ys[1]
end
νMax = 1.0
νs   = collect(range(0.0, νMax, dims[1])) # Positions of discretised vertices in polymerisation space 
dν   = νs[2]-νs[1]
nSpatialDims == 1 ? spacing  = [dν, dx] : spacing  = [dν, dx, dy]

A   = makeIncidenceMatrix3D(dims)
Ā   = abs.(A)
Aᵀ  = transpose(A)
Aᵤₚ = dropzeros((Ā-A).÷2)

# Number of edges over each dimension 
dimEdgeCount = Int64[]
for i=1:length(dims)
    push!(dimEdgeCount, (dims[i]-1)*prod(dims[Not(i)]))
end
nVerts  = prod(dims)          # Total number of vertices 
nEdges  = sum(dimEdgeCount)   # Total number of edges over all dimensions 

# Matrices for picking out ν and xy directions in derivatives 
Pν  = spdiagm(vcat(ones(Int64, dimEdgeCount[1]), zeros(Int64, sum(dimEdgeCount[2:end]))))   # Diagonal sparse matrix to exclude all xy edges 
Pxy  = spdiagm(vcat(zeros(Int64, dimEdgeCount[1]), ones(Int64, sum(dimEdgeCount[2:end]))))   # Diagonal sparse matrix to exclude all ν edges 

# Weights
W   = vertexVolumeWeightsMatrix(dims, spacing)
W⁻¹ =  vertexVolumeWeightsInverseMatrix(dims, spacing)
L⁻¹ = edgeLengthInverseMatrix(dims, spacing)

∇ₑ = L⁻¹*A       # Gradient operator giving gradient on each edge
∇cdot = -W⁻¹*Aᵀ  # Divergence operator giving divergence on each vertex calculated from edges 

# Diffusivity field over edges 
# Set no-flux boundary conditions by enforcing zero diffusivity in edges connection ghost points
Aperpₑ = edgePerpendicularAreaMatrix(dims, spacing)
𝓓ₑ     = 𝓓.*Aperpₑ # Sparse diagonal matrix of diffusivities over edges 

# Diagonal matrices of compartment thickness h over all vertices hᵥ
# Also diagonal matrix of thickness over edges, formed by taking mean of h at adjacent vertices 0.5.*Ā*hᵥ
thickness=="GRF" ? mat_h = hFun(dims, σ=0.1) : mat_h = ones(dims...)
hᵥ_vec = reshape(mat_h, nVerts)         # Cisternal thickness evaluated over vertices 
hₑ_vec = 0.5.*Ā*hᵥ_vec                  # Cisternal thickness evaluated over edges (mean of adjacent vertices)
hᵥ = spdiagm(hᵥ_vec)                    # Cisternal thickness over vertices, as a sparse diagonal matrix
hₑ = spdiagm(hₑ_vec)                    # Cisternal thickness over edges, as a sparse diagonal matrix
aᵥ = spdiagm(1.0./(1.0 .+ α_C.*hᵥ_vec)) # Prefactor 1/(1+α_C*hᵥ(x)) evaluated over vertices, packaged into a sparse diagonal matrix for convenience

uMat = zeros(Float64, dims...)
for ind in CartesianIndices(uMat)
    uMat[ind] = u0Fun([νs[ind[1]]], [0.0], [νMax/100.0])
end
# Ensure that the integral of concentration over ν at each point in space is 1
integ = spacing[1].*(0.5.*selectdim(uMat, 1, 1) .+ dropdims(sum(selectdim(uMat, 1, 2:dims[1]-1), dims=1), dims=1) .+ 0.5.*selectdim(uMat, 1, dims[1]))    

u0 = reshape(uMat, nVerts)
u0 .*= 1.0/integ[1]

# Set value of Fₑ at each point in space
# Integral of Fₑ over space is π
matFₑTmp = ones(Float64, dims[Not(1)]...)
for i=1:length(size(matFₑTmp))
    selectdim(matFₑTmp, i, 1) .*= 0.5
    selectdim(matFₑTmp, i, size(matFₑTmp)[i]) .*= 0.5
end
integF = prod(spacing[Not(1)])*sum(selectdim(matFₑTmp, 1, 1:size(matFₑTmp)[1]))

# Ensure integral of Fₑ over space is π
# matFₑ = (1/integF).*ones(Float64, dims[Not(1)]...)
matFₑ = ones(Float64, dims[Not(1)]...)
matE = zeros(dims...)
Esparse = spzeros(nVerts, nVerts)
E!(u0, dims, Esparse, matE, matFₑ, K₂, spacing[1])

# Part1 = aᵥ*∇cdot*Aperpₑ*(K₂*K₄.*Pν*∇ₑ - β.*Pν*Ā./2.0)
Part1 = aᵥ*∇cdot*Aperpₑ*(K₂*K₄.*Pν*∇ₑ - β.*Pν*Aᵤₚ)
Part2 = aᵥ*∇cdot*(hₑ*Pxy*𝓓ₑ*∇ₑ)

p = (Part1 = Part1, 
    Part2 = Part2, 
    u0 = u0, 
    dims = dims, 
    Esparse = Esparse, 
    matE = matE, 
    matFₑ = matFₑ, 
    K₂ = K₂, 
    dν = dν,
)
fullOperator = MatrixOperator(Esparse*Part1, update_func! = updateOperator!)
prob = ODEProblem(fullOperator, u0, (0.0, Tᵣ), p)

#%%

println("solving")
sol = solve(prob, Vern9(), saveat=Tᵣ/500.0, progress=true)

println("finished sim")

#%%
 
# Create directory for run data labelled with current time.
# paramsName = @savename nSpatialDims K₂ K₃ K₄ α_C δ_C σ N β 𝓓 Tᵣ h₀ Ωperp 𝓒
paramsName = @savename nSpatialDims K₂ K₄ α_C β 𝓓 Tᵣ
folderName = "$(Dates.format(Dates.now(),"yy-mm-dd-HH-MM-SS"))_$(paramsName)"
# Create frames subdirectory to store system state at each output time
subFolder = ""
mkpath(datadir("sims",subFolder,folderName))

#%%
spacing = [π^(1/nSpatialDims)/(Ngrid-1), 1/(Ngrid-1)]
W = vertexVolumeWeightsMatrix(dims, spacing)
if nSpatialDims==1
    concentrationSurfaceMovie(sol.u, sol.t, dims; subFolder=subFolder, folderName=folderName)
    spaceIntegralOver_ν_Movie(sol.u, sol.t, xs, νs, dims, W; subFolder=subFolder, folderName=folderName)
else
    spaceIntegralOver_ν_Movie(sol.u, sol.t, xs, νs, dims, W; subFolder=subFolder, folderName=folderName)
    uSlices = [selectdim(reshape(u, dims...), 3, dims[3]÷2) for u in sol.u]
    uSlicesReshaped = [reshape(u, prod(dims[Not(3)])) for u in uSlices]
    concentrationSurfaceMovie(uSlicesReshaped, sol.t, dims[1:2]; subFolder=subFolder, folderName=folderName)
end
