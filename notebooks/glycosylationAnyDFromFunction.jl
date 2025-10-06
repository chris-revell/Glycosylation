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
using JLD2

@from "$(srcdir("Glycosylation.jl"))" using Glycosylation
@from "$(srcdir("Visualise.jl"))" using Visualise
@from "$(srcdir("UsefulFunctions.jl"))" using UsefulFunctions
@from "$(srcdir("MakeWeightMatrices.jl"))" using MakeWeightMatrices
@from "$(srcdir("DerivedParameters.jl"))" using DerivedParameters

#%%

subFolder = "keynote"
terminateAt = "nuWall"
thicknessProfile = "GRF"
differencing = "centre"
solver = SSPRK432()
nOutputs = 100
σGRF = 0.3
λGRF = 0.05

nSpatialDims = 1
Ngrid = 201
dims = fill(Ngrid, nSpatialDims+1)

include(projectdir("notebooks", "paramsRaw.jl"))

#%%

derivedParams = derivedParameters(Ω, 𝒜, N, k_Cd, k_Ca, k_Sd, k_Sa, k₁, k₂, k₃, k₄, 𝒞, 𝒮, ℰ, D_C, D_S, Tᵣstar; checks=true)
@unpack L₀, E₀, C_b, S_b, δ_C, δ_S, α_C, α_S, C₀, S₀, Tᵣ, T̃ᵣ, K₂, K₃, K₄, σ, ϵ, 𝒟, β, h_C, h_S, λ, ζ, γ, Δ, F = derivedParams

#%%

# Create directory for run data labelled with current time.
paramsName = @savename nSpatialDims K₂ K₄ α_C β 𝒟 T̃ᵣ thicknessProfile differencing
folderName = "$(Dates.format(Dates.now(),"yy-mm-dd-HH-MM-SS"))_$(paramsName)"
# Create frames subdirectory to store system state at each output time
mkpath(datadir("sims",subFolder,folderName))

#%%

sol, p = glycosylationAnyD(dims, K₂, K₄, T̃ᵣ, α_C, 𝒟, β, thickness=thicknessProfile, differencing=differencing, solver=solver, nOutputs=nOutputs, σGRF=σGRF, λGRF=λGRF, terminateAt=terminateAt)
println("finished sim")

#%%

rawParams = (
    thicknessProfile = thicknessProfile,
    differencing = differencing,
    solver = solver,
    nOutputs = nOutputs,
    σGRF = σGRF,
    λGRF = λGRF,
    nSpatialDims = nSpatialDims,
    Ngrid = Ngrid,
    dims = dims,
    h₀ = h₀,
    𝒜 = 𝒜,
    Ω = Ω,
    N = N,
    k_Cd = k_Cd,
    k_Ca = k_Ca,
    k_Sd = k_Sd,
    k_Sa = k_Sa,
    k₁ = k₁,
    k₂ = k₂,
    k₃ = k₃,
    k₄ = k₄,
    𝒞 = 𝒞,
    𝒮 = 𝒮,
    ℰ = ℰ,
    D_C = D_C,
    D_S = D_S,
    Tᵣstar = Tᵣstar,
    ϕ = ϕ
)

jldsave(datadir("sims",subFolder,folderName,"solution.jld2"); sol, p, rawParams)

#%%

if nSpatialDims==1
    concentrationSurfaceMovie(sol.u, dims; subFolder=subFolder, folderName=folderName)
    # concentrationHeatmapMovie(sol.u, dims; subFolder=subFolder, folderName=folderName)
    M̃movie(sol.u, p; subFolder=subFolder, folderName=folderName)
    if thicknessProfile=="GRF"
        thicknessPlot(p.hᵥ, p.dims; subFolder=subFolder, folderName=folderName)
    end
else    
    uSlices = [selectdim(reshape(u, dims...), 3, dims[3]÷2) for u in sol.u[1:35]]
    uSlicesReshaped = [reshape(u, prod(dims[Not(3)])) for u in uSlices[1:35]]
    concentrationSurfaceMovie(uSlicesReshaped, dims[1:2]; subFolder=subFolder, folderName=folderName)
    # concentrationHeatmapMovie(uSlicesReshaped, dims; subFolder=subFolder, folderName=folderName)
    M̃movie(sol.u[1:50], p; subFolder=subFolder, folderName=folderName)
    if thicknessProfile=="GRF"
        thicknessPlot(p.hᵥ, dims; subFolder=subFolder, folderName=folderName)
    end
end
