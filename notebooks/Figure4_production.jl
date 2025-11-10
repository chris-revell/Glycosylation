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
@from "$(srcdir("UsefulFunctions.jl"))" using UsefulFunctions

include(projectdir("notebooks", "paramsRaw.jl"))

subFolder = "Figure4"
terminateAt = "halfProduction"
σGaussian = 0.20
nSpatialDims = 1
Ngrid = 401
dims = fill(Ngrid, nSpatialDims+1)

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
folderName = "$(Dates.format(Dates.now(),"yy-mm-dd-HH-MM-SS"))_$(paramsName)_production"
mkpath(datadir("sims",subFolder,folderName))

#%%
sol0, p0 = glycosylation(dims, K₂, K₄, T̃ᵣ, α_C, 𝒟, β, thickness="uniform", differencing=differencing, solver=solver, nOutputs=nOutputs, terminateAt=terminateAt, saveIntermediate=false)
println("finished sim 0")
jldsave(datadir("sims",subFolder,folderName,"solutionNoVariation.jld2"); sol0, p0)#, rawParams1)

Tᵣ₅₀Star = sol0.t[end]*(N^2)*(K₂+σ*K₃)/(k₁*E₀)
𝒫simNoVariation = Mstarϕ(sol0.u[end], p0.W, p0.dims, p0.dν, p0.hᵥ, α_C, 𝒞, 0.5)/Tᵣ₅₀Star
@show 𝒫simNoVariation


sol1, p1 = glycosylation(dims, K₂, K₄, T̃ᵣ, α_C, 𝒟, β, thickness="Gaussian", differencing=differencing, solver=solver, nOutputs=nOutputs, σGaussian=σGaussian, terminateAt=terminateAt, saveIntermediate=false)
println("finished sim 1")
jldsave(datadir("sims",subFolder,folderName,"solutionHVariation.jld2"); sol1, p1)#, rawParams1)

Tᵣ₅₀Star = sol1.t[end]*(N^2)*(K₂+σ*K₃)/(k₁*E₀)
𝒫simSpatialVariation = Mstarϕ(sol1.u[end], p1.W, p1.dims, p1.dν, p1.hᵥ, α_C, 𝒞, 0.5)/Tᵣ₅₀Star
@show 𝒫simSpatialVariation

#%%

sol2, p2 = glycosylation(dims, K₂, K₄, T̃ᵣ, α_C, 𝒟, β, thickness="uniform", fDist="Gaussian", differencing=differencing, solver=solver, nOutputs=nOutputs, σGaussian=σGaussian, terminateAt=terminateAt, saveIntermediate=false)
println("finished sim 2")
jldsave(datadir("sims",subFolder,folderName,"solutionFVariation.jld2"); sol2, p2)#, rawParams2)

Tᵣ₅₀Star = sol2.t[end]*(N^2)*(K₂+σ*K₃)/(k₁*E₀)
𝒫simEnzymeVariation = Mstarϕ(sol2.u[end], p2.W, p2.dims, p2.dν, p2.hᵥ, α_C, 𝒞, 0.5)/Tᵣ₅₀Star
@show 𝒫simEnzymeVariation