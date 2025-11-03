using OrdinaryDiffEq
using SparseArrays
using UnPack
using FromFile
using DrWatson
using Printf
using SciMLOperators
using Dates
using InvertedIndices
using Statistics
using JLD2

@from "$(srcdir("Glycosylation.jl"))" using Glycosylation
@from "$(srcdir("UsefulFunctions.jl"))" using UsefulFunctions
@from "$(srcdir("MakeWeightMatrices.jl"))" using MakeWeightMatrices
@from "$(srcdir("DerivedParameters.jl"))" using DerivedParameters

include(projectdir("notebooks", "paramsRaw.jl"))
terminateAt = "nuWall"
thicknessProfile = "uniform"
nSpatialDims = 1
Ngrid = 401
dims = fill(Ngrid, nSpatialDims+1)

subFolder = "Figure3"
folderName = "$(Dates.format(Dates.now(),"yy-mm-dd-HH-MM-SS"))" # Create directory for run data labelled with current time.
mkpath(datadir("sims",subFolder,folderName))

h_C = 2*k_Ca/k_Cd
h_S = 2*k_Sa/k_Sd
hcutoff = (2.0*k_Sa/k_Sd)*((𝒮*k₁*k₃)/(2.0*𝒜*k₂*k₄) - 1.0)
h₀s = collect(h_C/10.0:h_C/2.0:h_C*5)
# h₀s2 = collect(hMax+1.0:1.0:7.0)
h₀s2 = collect(h₀s[end]+0.5:h_C*2:hcutoff-0.001)
append!(h₀s, h₀s2)

Ωs = h₀s.*𝒜      # Dimensional lumen volume 

sols = []
ps = []
for i=1:length(h₀s)
    @show h₀s[i]    
    derivedParams = derivedParameters(Ωs[i], 𝒜, N, k_Cd, k_Ca, k_Sd, k_Sa, k₁, k₂, k₃, k₄, 𝒞, 𝒮, ℰ, D_C, D_S, Tᵣstar; checks=false)
    @unpack L₀, E₀, C_b, S_b, δ_C, δ_S, α_C, α_S, C₀, S₀, Tᵣ, T̃ᵣ, K₂, K₃, K₄, σ, ϵ, 𝒟, β, h_C, h_S, λ, ζ, γ, Δ, F = derivedParams
    sol, p = glycosylation(dims, K₂, K₄, 10.0*T̃ᵣ, α_C, 𝒟, β, thickness=thicknessProfile, differencing=differencing, solver=solver, nOutputs=nOutputs, terminateAt="halfProduction", saveIntermediate=false) 
    Tᵣ₅₀Star = sol.t[end]*(N^2)*(K₂+σ*K₃)/(k₁*E₀)
    push!(sols, sol)
    push!(ps, p)
end

jldsave(datadir("sims",subFolder,folderName,"solutions.jld2"); sols, ps, h₀s)
