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
using Statistics

@from "$(srcdir("Glycosylation.jl"))" using Glycosylation
@from "$(srcdir("Visualise.jl"))" using Visualise
@from "$(srcdir("UsefulFunctions.jl"))" using UsefulFunctions
@from "$(srcdir("MakeWeightMatrices.jl"))" using MakeWeightMatrices
@from "$(srcdir("DerivedParameters.jl"))" using DerivedParameters
@from "$(srcdir("MakeIncidenceMatrix.jl"))" using MakeIncidenceMatrix

u0fun(xs, μs, σs) = exp(-sum((xs.-μs).^2.0./σs.^2.0)) # Multidimensional Gaussian


nSpatialDims = 1

# h₀s = collect(0.000001:0.0000005:0.00001)
# h₀ = h₀s[1]
h₀ = 1.0

Ωperp = 1.0  # Lumen footprint area
N     = 100         # Maximum polymer length 
# k_Cd  = 300000.0 # Complex desorption rate
# k_Ca  = 0.1 # Complex adsorption rate
# k_Sd  = 1.0 # Substrate desorption rate
# k_Sa  = 1.0 # Substrate adsorption rate
# k₁    = 2.0   # Complex formation forward reaction rate 
# # k₂    = 0.01   # Complex dissociation reverse reaction rate 
# k₂    = 1.0   # Complex dissociation reverse reaction rate 
# k₃    = 0.003   # Product formation
# k₄    = 2.0  # Product dissociation 
# E_0   = 0.01
# 𝓒     = 1.0
# 𝓢     = 1000.0
# D_C   = 0.00000001  # Monomer/polymer diffusivity
# D_S   = 0.00000001  # Substrate diffusivity
# Tᵣstar= 0.1  # Release time
# ϕ     = 0.5

Ngrid = 101
nSpatialDims == 1 ? dims  = [Ngrid, Ngrid] : dims  = [Ngrid, Ngrid, Ngrid]
# derivedParams = derivedParameters(h₀, Ωperp, N, k_Cd, k_Ca, k_Sd, k_Sa, k₁, k₂, k₃, k₄, E_0, 𝓒, 𝓢, D_C, D_S, Tᵣstar; checks=true)
# @unpack 𝓔, K₃, K₄, δ_C, δ_S, Tᵣ, Ω, α_C, α_S, C_b, S_b, C_0, S_0, K₂, σ, ϵ, 𝓓, β, K₂, L₀ = derivedParams

𝓒 = 1.0
K₂ = 1.0
K₄ = 1.0
Tᵣ = 0.5
α_C = 1.0
𝓓 = 1.0
β = 1.0

#%%

xMax = (Ωperp)^(1/nSpatialDims)
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

mat_h = h₀.*ones(fill(Ngrid, nSpatialDims+1)...)

sol = glycosylationAnyD(mat_h, dims, Ωperp, 𝓒, K₂, K₄, Tᵣ, α_C , 𝓓, β)


# # PDE discretisation parameters 
# nSpatialDims = length(dims)-1
    
# xMax = (Ωperp)^(1/nSpatialDims)
# xs   = collect(range(0.0, xMax, dims[2]))
# dx   = xs[2]-xs[1]
# if nSpatialDims > 1 
#     yMax = xMax
#     ys   = collect(range(0.0, yMax, dims[3]))
#     dy   = ys[2]-ys[1]
# end
# νMax = 1.0
# νs   = collect(range(0.0, νMax, dims[1])) # Positions of discretised vertices in polymerisation space 
# dν   = νs[2]-νs[1]
# nSpatialDims == 1 ? spacing  = [dν, dx] : spacing  = [dν, dx, dy]

# h₀ = mean(selectdim(mat_h, 1, 1))

# A   = makeIncidenceMatrix3D(dims)
# Ā   = abs.(A)
# Aᵀ  = transpose(A)
# Aᵤₚ = dropzeros((Ā-A).÷2)

# # Number of edges over each dimension 
# dimEdgeCount = Int64[]
# for i=1:length(dims)
#     push!(dimEdgeCount, (dims[i]-1)*prod(dims[Not(i)]))
# end
# nVerts  = prod(dims)          # Total number of vertices 
# nEdges  = sum(dimEdgeCount)   # Total number of edges over all dimensions 

# # Matrices for picking out ν and xy directions in derivatives 
# Pν  = spdiagm(vcat(ones(Int64, dimEdgeCount[1]), zeros(Int64, sum(dimEdgeCount[2:end]))))   # Diagonal sparse matrix to exclude all xy edges 
# Pxy  = spdiagm(vcat(zeros(Int64, dimEdgeCount[1]), ones(Int64, sum(dimEdgeCount[2:end]))))   # Diagonal sparse matrix to exclude all ν edges 

# # Weights
# W   = vertexVolumeWeightsMatrix(dims, spacing)
# W⁻¹ =  vertexVolumeWeightsInverseMatrix(dims, spacing)
# L⁻¹ = edgeLengthInverseMatrix(dims, spacing)

# ∇ₑ = L⁻¹*A       # Gradient operator giving gradient on each edge
# ∇cdot = -W⁻¹*Aᵀ  # Divergence operator giving divergence on each vertex calculated from edges 

# # Diffusivity field over edges 
# # Set no-flux boundary conditions by enforcing zero diffusivity in edges connection ghost points
# Aperpₑ = edgePerpendicularAreaMatrix(dims, spacing)
# 𝓓ₑ     = 𝓓.*Aperpₑ # Sparse diagonal matrix of diffusivities over edges 

# # Diagonal matrices of compartment thickness h over all vertices hᵥ
# # Also diagonal matrix of thickness over edges, formed by taking mean of h at adjacent vertices 0.5.*Ā*hᵥ
# hᵥ_vec = reshape(mat_h, nVerts)         # Cisternal thickness evaluated over vertices 
# hₑ_vec = 0.5.*Ā*hᵥ_vec                  # Cisternal thickness evaluated over edges (mean of adjacent vertices)
# hᵥ = spdiagm(hᵥ_vec)                    # Cisternal thickness over vertices, as a sparse diagonal matrix
# hₑ = spdiagm(hₑ_vec)                    # Cisternal thickness over edges, as a sparse diagonal matrix
# aᵥ = spdiagm(1.0./(1.0 .+ α_C.*hᵥ_vec)) # Prefactor 1/(1+α_C*hᵥ(x)) evaluated over vertices, packaged into a sparse diagonal matrix for convenience

# uMat = zeros(Float64, dims...)
# for ind in CartesianIndices(uMat)
#     uMat[ind] = u0fun([νs[ind[1]]], [0.0], [νMax/100.0])
# end
# integ = 0.5*spacing[1]*(uMat[1,1]+uMat[end,1]) + spacing[1]*sum(uMat[2:end-1,1])
# u0 = reshape(uMat, nVerts)
# # integ = sum(W*hᵥ*u0)
# u0 .*= 1/integ

# # Set value of Fₑ at each point in space
# matFₑTmp = ones(Float64, dims[Not(1)]...)
# for i=1:length(size(matFₑTmp))
#     selectdim(matFₑTmp, i, 1) .*= 0.5
#     selectdim(matFₑTmp, i, size(matFₑTmp)[i]) .*= 0.5
# end
# integF = prod(spacing[Not(1)])*sum(selectdim(matFₑTmp, 1, 1:size(matFₑTmp)[1]))
# # Ensure integral of Fₑ over space is π
# matFₑ = (π/integF).*ones(Float64, dims[Not(1)]...)

# matE = zeros(dims...)
# Esparse = spzeros(nVerts, nVerts)
# E!(u0, dims, Esparse, matE, matFₑ, K₂, dν)

# # PDE operator components
# # Part1 = aᵥ*∇cdot*(K₂*K₄.*Pν*∇ₑ - β.*Pν*Aᵤₚ)
# # Part2 = 𝓓.*aᵥ*∇cdot*(hₑ*Pxy*∇ₑ)

# Part1 = aᵥ*∇cdot*Aperpₑ*(K₂*K₄.*Pν*∇ₑ - β.*Pν*Aᵤₚ)
# Part2 = aᵥ*∇cdot*(hₑ*Pxy*𝓓ₑ*∇ₑ)

# p = (Part1 = Part1, 
#     Part2 = Part2, 
#     u0 = u0, 
#     dims = dims, 
#     Esparse = Esparse, 
#     matE = matE, 
#     matFₑ = matFₑ, 
#     K₂ = K₂, 
#     dν = dν,
# )
# fullOperator = MatrixOperator(Esparse*Part1.+Part2, update_func! = updateOperator!)
# prob = ODEProblem(fullOperator, u0, (0.0, Tᵣ), p)


# #%%


# println("solving")




# sol = solve(prob, Vern9(), saveat=Tᵣ/50.0, progress=true)






# println("finished sim")

#%%

# Create directory for run data labelled with current time.
paramsName = @savename nSpatialDims K₂ K₄ α_C N β 𝓓 Tᵣ h₀ Ωperp 𝓒
folderName = "$(Dates.format(Dates.now(),"yy-mm-dd-HH-MM-SS"))_$(paramsName)"
# Create frames subdirectory to store system state at each output time
subFolder = ""
mkpath(datadir("sims",subFolder,folderName))

W = vertexVolumeWeightsMatrix(dims, spacing)

# if nSpatialDims==1
#     concentrationSurfaceMovie(sol.u, sol.t, dims; subFolder=subFolder, folderName=folderName)
#     spaceIntegralOver_ν_Movie(sol.u, sol.t, xs, νs, dims, W; subFolder=subFolder, folderName=folderName)
# else
#     spaceIntegralOver_ν_Movie(sol.u, sol.t, xs, νs, dims, W; subFolder=subFolder, folderName=folderName)
#     uSlices = [reshape(u, dims...)[:,:,dims[3]÷2] for u in sol.u]
#     uSlicesReshaped = [reshape(u, prod(dims[Not(3)])) for u in uSlices]
#     concentrationSurfaceMovie(uSlicesReshaped, sol.t, dims[1:2]; subFolder=subFolder, folderName=folderName)
# end

# Ẽ = π*K₂/(1+K₂)
# a = Ẽ*β*Tᵣ
# b = 1+α_C

# println("$(a < b)")

#%%

fig = Figure(size=(1000,1000))
ax = CairoMakie.Axis(fig[1, 1], aspect=1)
ax.xlabel = "ν"
ax.ylabel = "C"
analyticLine = Observable(zeros(dims[1]))
numericLine = Observable(zeros(dims[1]))
lines!(ax, νs, analyticLine, color=:red)
lines!(ax, νs, numericLine, color=:blue)
# maxY = [0.0]
ylims!(ax, (0.0, 15.0))
record(fig, datadir("sims",subFolder, folderName, "analyticCs.mp4"), 1:length(sol.t[2:end]); framerate=10) do i
    tst = homogeneousWidthC.(νs, K₂, K₄, α_C, β, sol.t[i])
    # @show maximum(tst)
    # maxY[1] = maximum([maxY[1], maximum(tst)])
    analyticLine[] .= tst
    uInternal = reshape(sol.u[i], dims...)
    numericLine[] .= 10.0.*uInternal[:,dims[2]÷2]
    analyticLine[] = analyticLine[]
    numericLine[] = numericLine[]
    
end


