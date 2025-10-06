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

subFolder = "Figure2"
# Create directory for run data labelled with current time.
folderName = "$(Dates.format(Dates.now(),"yy-mm-dd-HH-MM-SS"))"
# Create frames subdirectory to store system state at each output time
mkpath(datadir("sims",subFolder,folderName))
terminateAt = "nuWall"
thicknessProfile = "uniform"
differencing = "centre"
solver = SSPRK432()
nOutputs = 10000
σGRF = 0.2

nSpatialDims = 1
Ngrid = 201
dims = fill(Ngrid, nSpatialDims+1)

include(projectdir("notebooks", "paramsRaw.jl"))

#%%

h_C = 2*k_Ca/k_Cd
h_S = 2*k_Sa/k_Sd
hMax = h_C*5
hMin = h_C/10
h₀s = collect(hMin:2*hMin:hMax)

h₀s2 = collect(hMax+1.0:1.0:7.0)
append!(h₀s, h₀s2)

Ωs = h₀s.*𝒜      # Dimensional lumen volume 

#%%

sols = []
ps = []
for i=1:length(h₀s)
    @show h₀s[i]    
    derivedParams = derivedParameters(Ωs[i], 𝒜, N, k_Cd, k_Ca, k_Sd, k_Sa, k₁, k₂, k₃, k₄, 𝒞, 𝒮, ℰ, D_C, D_S, Tᵣstar; checks=false)
    @unpack L₀, E₀, C_b, S_b, δ_C, δ_S, α_C, α_S, C₀, S₀, Tᵣ, T̃ᵣ, K₂, K₃, K₄, σ, ϵ, 𝒟, β, h_C, h_S, λ, ζ, γ, Δ, F = derivedParams
    sol, p = glycosylationAnyD(dims, K₂, K₄, 1000.0, α_C, 𝒟, β, thickness=thicknessProfile, differencing=differencing, solver=solver, nOutputs=nOutputs, terminateAt="halfProduction", saveIntermediate=false) 
    Tᵣ₅₀Star = sol.t[end]*(N^2)*(K₂+σ*K₃)/(k₁*E₀)
    push!(sols, sol)
    push!(ps, p)
end

#%%

derivedParams = derivedParameters(Ω, 𝒜, N, k_Cd, k_Ca, k_Sd, k_Sa, k₁, k₂, k₃, k₄, 𝒞, 𝒮, ℰ, D_C, D_S, Tᵣstar; checks=false)
@unpack L₀, E₀, C_b, S_b, δ_C, δ_S, α_C, α_S, C₀, S₀, Tᵣ, T̃ᵣ, K₂, K₃, K₄, σ, ϵ, 𝒟, β, h_C, h_S, λ, ζ, γ, Δ, F = derivedParams
sol1, p1 = glycosylationAnyD(dims, K₂, K₄, T̃ᵣ, α_C, 𝒟, β, thickness=thicknessProfile, differencing=differencing, solver=solver, nOutputs=nOutputs, terminateAt=terminateAt)
println("finished sim")

#%%

