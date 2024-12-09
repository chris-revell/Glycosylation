
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

@from "$(srcdir("Glycosylation.jl"))" using Glycosylation
@from "$(srcdir("Visualise.jl"))" using Visualise
@from "$(srcdir("UsefulFunctions.jl"))" using UsefulFunctions
@from "$(srcdir("MakeWeightMatrices.jl"))" using MakeWeightMatrices
@from "$(srcdir("DerivedParameters.jl"))" using DerivedParameters
@from "$(srcdir("CisternaWidth.jl"))" using CisternaWidth

#%%

thicknessProfile = "uniform"
differencing = "centre"
solver = SSPRK432()
nOutputs = 100
σ = 0.1


nSpatialDims = 2
Ngrid = 201
# dims = [Ngrid,2]
dims = fill(Ngrid, nSpatialDims+1)
# dims[1] = 01
#%%

h₀ = 0.1
Ωperp = 10000    # Dimensional lumen footprint area
Ω     = h₀*Ωperp      # Dimensional lumen volume 
N     = 100     # Maximum polymer length 
k_Cd  = 1.0 # Complex desorption rate
k_Ca  = 0.01 # Complex adsorption rate
k_Sd  = 1.0 # Substrate desorption rate
k_Sa  = 0.01 # Substrate adsorption rate
k₁    = 1.0   # Complex formation forward reaction rate 
k₂    = 0.1   # Complex dissociation reverse reaction rate 
k₃    = 0.1   # Product formation
k₄    = 0.1  # Product dissociation 
𝓒     = 100000.0
𝓢     = 100000.0
𝓔     = 0.0001
D_C   = 0.0000001  # Monomer/polymer diffusivity
D_S   = 0.0000001  # Substrate diffusivity
Tᵣstar= 1000000000.0  # Release time
ϕ     = 0.5

#%%

derivedParams = derivedParameters(Ω, Ωperp, N, k_Cd, k_Ca, k_Sd, k_Sa, k₁, k₂, k₃, k₄, 𝓒, 𝓢, 𝓔, D_C, D_S, Tᵣstar; checks=true)
@unpack L₀, E₀, h₀, C_b, S_b, δ_C, δ_S, α_C, α_S, C₀, S₀, Tᵣ, K₂, K₃, K₄, σ, ϵ, 𝓓, β, λ = derivedParams

#%%

sol, mat_h = glycosylationAnyD(dims, K₂, K₄, Tᵣ, α_C, 𝓓, β, thickness=thicknessProfile, differencing=differencing, solver=solver, nOutputs=nOutputs, σ=σ)
println("finished sim")

#%%

# Create directory for run data labelled with current time.
# paramsName = @savename nSpatialDims K₂ K₃ K₄ α_C δ_C σ N β 𝓓 Tᵣ h₀ Ωperp 𝓒
paramsName = @savename nSpatialDims K₂ K₄ α_C β 𝓓 Tᵣ thicknessProfile differencing
folderName = "$(Dates.format(Dates.now(),"yy-mm-dd-HH-MM-SS"))_$(paramsName)"
# Create frames subdirectory to store system state at each output time
subFolder = ""
mkpath(datadir("sims",subFolder,folderName))

#%%

xMax = 0.01*π^(1/nSpatialDims)
xs   = collect(range(0.0, xMax, dims[2]))
dx   = xs[2]-xs[1]
if nSpatialDims > 1 
    yMax = xMax
    ys   = collect(range(0.0, yMax, dims[3]))
    dy   = ys[2]-ys[1]
end
νMax = 1.0
νs   = collect(range(0.0, νMax, dims[1]))
dν   = νs[2]-νs[1]
nSpatialDims == 1 ? spacing  = [dν, dx] : spacing  = [dν, dx, dy]
W = vertexVolumeWeightsMatrix(dims, spacing)

#%%

if nSpatialDims==1
    concentrationSurfaceMovie(sol.u, sol.t, dims; subFolder=subFolder, folderName=folderName)
    concentrationHeatmapMovie(sol.u, sol.t, dims; subFolder=subFolder, folderName=folderName)
    spaceIntegralOver_ν_Movie(sol.u, sol.t, xs, νs, dims, W; subFolder=subFolder, folderName=folderName)
    if thicknessProfile=="GRF"
        thicknessPlot(mat_h; subFolder=subFolder, folderName=folderName)
    end
else
    spaceIntegralOver_ν_Movie(sol.u, sol.t, xs, νs, dims, W; subFolder=subFolder, folderName=folderName)
    uSlices = [selectdim(reshape(u, dims...), 3, dims[3]÷2) for u in sol.u]
    uSlicesReshaped = [reshape(u, prod(dims[Not(3)])) for u in uSlices]
    concentrationSurfaceMovie(uSlicesReshaped, sol.t, dims[1:2]; subFolder=subFolder, folderName=folderName)
    concentrationHeatmapMovie(uSlicesReshaped, sol.t, dims[1:2]; subFolder=subFolder, folderName=folderName)
    if thicknessProfile=="GRF"
        thicknessPlot(mat_h; subFolder=subFolder, folderName=folderName)
    end
end


# W = vertexVolumeWeightsMatrix(dims, spacing)
# hᵥ = spdiagm(reshape(mat_h, prod(dims)))
# h₀ = 0.01
# Ωperp = 10000
# Ω = Ωperp*h₀
# N = 100
# 𝓒     = 100000.0
# 𝓢     = 100000.0
# 𝓔     = 0.0001
# C_b  = 𝓒/Ω
# k₁ = 1.0
# @show P_star(sol.u[end], W, dims, dν, hᵥ, α_C, C_b, Ω, ϕ, Ωperp, k₁, 𝓔, Tᵣ)
# @show 𝓟starUniform(𝓒, 𝓔, 𝓢, ϕ, N, k₁, k₂, K₃, K₄, Ωperp, h₀, h_C, h_S)

# Mϕ50 = α_C*𝓒/(2.0*(1+α_C))


# T50 = 2.0*Ωperp/(k₁*𝓔) * N^2* (K₂+σ*K₃) * (t_0 + (ϕ-ν_0)*(1+α_C)*())
# P50 = Mϕ50/T50


# Tᵣ = 30.0
# K₂ = 1.0
# K₄ = 0.0001
# α_C = 1.0
# 𝓓 = 1.0
# β = 0.1

# differencing = "upstream"

# include(projectdir("notebooks","paramsRaw.jl"))
# h_C = 2*k_Ca/k_Cd
# h_S = 2*k_Sa/k_Sd
# hMax = h_C*25
# hMin = h_C/10
# h₀s = collect(hMin:2*hMin:hMax)
# Ωs    = h₀s.*Ωperp      # Dimensional lumen volume 

# h₀ = 2*h_C
# Ω = h₀*Ωperp
