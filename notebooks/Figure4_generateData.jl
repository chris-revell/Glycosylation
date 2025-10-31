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
using JLD2
using LinearAlgebra

@from "$(srcdir("Glycosylation.jl"))" using Glycosylation
@from "$(srcdir("DerivedParameters.jl"))" using DerivedParameters

include(projectdir("notebooks", "paramsRaw.jl"))

subFolder = "new/Figure4"
terminateAt = "nuWall"
σGaussian = 0.20
nSpatialDims = 1
Ngrid = 201
dims = fill(Ngrid, nSpatialDims+1)

# K₂ = 0.3
# K₄ = 1.0
# T̃ᵣ = 0.4
# α_C = 5.0
# 𝒟 = 200.0
# β = 70.0

# h₀ = 1.0

rawParams1 = (
    thicknessProfile = "Gaussian",
    differencing = differencing,
    solver = solver,
    nOutputs = nOutputs,
    σGaussian = σGaussian,
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
rawParams2 = (
    thicknessProfile = "uniform",
    differencing = differencing,
    solver = solver,
    nOutputs = nOutputs,
    σGaussian = σGaussian,
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

derivedParams = derivedParameters(Ω, 𝒜, N, k_Cd, k_Ca, k_Sd, k_Sa, k₁, k₂, k₃, k₄, 𝒞, 𝒮, ℰ, D_C, D_S, Tᵣstar; checks=true)
@unpack L₀, E₀, C_b, S_b, δ_C, δ_S, α_C, α_S, C₀, S₀, Tᵣ, T̃ᵣ, K₂, K₃, K₄, σ, ϵ, 𝒟, β, h_C, h_S, λ, ζ, γ, Δ, F = derivedParams

#%%

# Create directory for run data labelled with current time.
paramsName = @savename nSpatialDims K₂ K₄ α_C β 𝒟 T̃ᵣ differencing
folderName = "$(Dates.format(Dates.now(),"yy-mm-dd-HH-MM-SS"))_$(paramsName)"
mkpath(datadir("sims",subFolder,folderName))

#%%

sol1, p1 = glycosylation(dims, K₂, K₄, T̃ᵣ, α_C, 𝒟, β, thickness="Gaussian", differencing=differencing, solver=solver, nOutputs=nOutputs, σGaussian=σGaussian, terminateAt=terminateAt)
println("finished sim 1")
jldsave(datadir("sims",subFolder,folderName,"solutionHVariation.jld2"); sol1, p1)#, rawParams1)

#%%

sol2, p2 = glycosylation(dims, K₂, K₄, T̃ᵣ, α_C, 𝒟, β, thickness="uniform", fDist="Gaussian", differencing=differencing, solver=solver, nOutputs=nOutputs, σGaussian=σGaussian, terminateAt=terminateAt)
println("finished sim 2")
jldsave(datadir("sims",subFolder,folderName,"solutionFVariation.jld2"); sol2, p2)#, rawParams2)

