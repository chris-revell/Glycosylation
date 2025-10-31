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
# using XLSX
# using DataFrames
# using Interpolations
using Statistics
using JLD2

@from "$(srcdir("Glycosylation.jl"))" using Glycosylation
@from "$(srcdir("Visualise.jl"))" using Visualise
@from "$(srcdir("UsefulFunctions.jl"))" using UsefulFunctions
@from "$(srcdir("MakeWeightMatrices.jl"))" using MakeWeightMatrices
@from "$(srcdir("DerivedParameters.jl"))" using DerivedParameters

#%%

subFolder = ""
terminateAt = "nuWall"
thicknessProfile = "GRF"
# nOutputs = 100
σGRF = 0.3
λGRF = 0.1

nSpatialDims = 2
Ngrid = 101
dims = fill(Ngrid, nSpatialDims+1)

params = "raw"

if params=="raw"
    include(projectdir("notebooks", "paramsRaw.jl"))
    derivedParams = derivedParameters(Ω, 𝒜, N, k_Cd, k_Ca, k_Sd, k_Sa, k₁, k₂, k₃, k₄, 𝒞, 𝒮, ℰ, D_C, D_S, Tᵣstar; checks=true)
    @unpack L₀, E₀, C_b, S_b, δ_C, δ_S, α_C, α_S, C₀, S₀, Tᵣ, T̃ᵣ, K₂, K₃, K₄, σ, ϵ, 𝒟, β, h_C, h_S, λ, ζ, γ, Δ, F = derivedParams
    savedParams = (
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
else
    include(projectdir("notebooks", "paramsDerived.jl"))
    savedParams = (
        K₂ = K₂,
        K₄ = K₄,
        T̃ᵣ = T̃ᵣ,
        α_C = α_C,
        𝒟 = 𝒟,
        β = β,
    )
end
#%%

# Create directory for run data labelled with current time.
paramsName = @savename nSpatialDims K₂ K₄ α_C β 𝒟 T̃ᵣ thicknessProfile differencing
folderName = "$(Dates.format(Dates.now(),"yy-mm-dd-HH-MM-SS"))_$(paramsName)"
mkpath(datadir("sims",subFolder,folderName))

#%%

sol, p = glycosylation(dims, K₂, K₄, T̃ᵣ, α_C, 𝒟, β, thickness=thicknessProfile, differencing=differencing, solver=solver, nOutputs=nOutputs, σGRF=σGRF, λGRF=λGRF, terminateAt=terminateAt)
println("finished sim")

jldsave(datadir("sims",subFolder,folderName,"solution.jld2"); sol, p, savedParams)

#%%

if nSpatialDims==1
    concentrationSurfaceMovie(sol.u, dims; subFolder=subFolder, folderName=folderName)
    # concentrationHeatmapMovie(sol.u, dims; subFolder=subFolder, folderName=folderName)
    M̃movie(sol.u, p; subFolder=subFolder, folderName=folderName)
    if thicknessProfile=="GRF"
        thicknessPlot(p.hᵥ, p.dims; subFolder=subFolder, folderName=folderName)
    end
else    
    uSlices = [selectdim(reshape(u, dims...), 3, dims[3]÷2) for u in sol.u]
    uSlicesReshaped = [reshape(u, prod(dims[Not(3)])) for u in uSlices]
    concentrationSurfaceMovie(uSlicesReshaped, dims[1:2]; subFolder=subFolder, folderName=folderName)
    # concentrationHeatmapMovie(uSlicesReshaped, dims; subFolder=subFolder, folderName=folderName)
    M̃movie(sol.u, p; subFolder=subFolder, folderName=folderName)
    if thicknessProfile=="GRF"
        thicknessPlot(p.hᵥ, dims; subFolder=subFolder, folderName=folderName)
    end
end
