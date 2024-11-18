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

@from "$(srcdir("Glycosylation.jl"))" using Glycosylation
@from "$(srcdir("Visualise.jl"))" using Visualise
@from "$(srcdir("UsefulFunctions.jl"))" using UsefulFunctions
@from "$(srcdir("MakeWeightMatrices.jl"))" using MakeWeightMatrices
@from "$(srcdir("DerivedParameters.jl"))" using DerivedParameters
@from "$(srcdir("CisternaWidth.jl"))" using CisternaWidth

#%%

thicknessProfile = "GRF"

include(projectdir("notebooks","paramsRaw.jl"))
h_C = 2*k_Ca/k_Cd
h_S = 2*k_Sa/k_Sd
hMax = h_C*25
hMin = h_C/10
h₀s = collect(hMin:2*hMin:hMax)
Ωs    = h₀s.*Ωperp      # Dimensional lumen volume 

h₀ = 2*h_C
Ω = h₀*Ωperp

nSpatialDims = 2
Ngrid = 51
dims  = fill(Ngrid, nSpatialDims+1)

derivedParams = derivedParameters(Ω, Ωperp, N, k_Cd, k_Ca, k_Sd, k_Sa, k₁, k₂, k₃, k₄, 𝓒, 𝓢, 𝓔, D_C, D_S, Tᵣstar; checks=true)
@unpack L₀, E₀, h₀, C_b, S_b, δ_C, δ_S, α_C, α_S, C₀, S₀, Tᵣ, K₂, K₃, K₄, σ, ϵ, 𝓓, β, λ = derivedParams

#%%

# include(projectdir("notebooks", "paramsDerived.jl"))

#%%

sol, mat_h = glycosylationAnyD(dims, K₂, K₄, Tᵣ, α_C, 𝓓, β, thickness=thicknessProfile, differencing="upstream", solver=SSPRK432(), nOutputs=50)
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
xMax = π^(1/nSpatialDims)
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