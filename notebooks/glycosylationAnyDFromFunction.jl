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

nSpatialDims = 1
Ngrid = 201

#%%

Ωperp = 1000.0    # Dimensional lumen footprint area
Ω     = 1.0      # Dimensional lumen volume 
N     = 1000     # Maximum polymer length 
k_Cd  = 100.0    # Dimensional complex desorption rate
k_Ca  = 0.01     # Dimensional complex adsorption rate
k_Sd  = 1.0      # Dimensional substrate desorption rate
k_Sa  = 1.0      # Dimensional substrate adsorption rate
k₁    = 1.0      # Dimensional complex formation forward reaction rate 
k₂    = 1.0     # Dimensional complex dissociation reverse reaction rate 
k₃    = 1.0     # Dimensional product formation
k₄    = 1.0      # Dimensional product dissociation 
𝓔     = 0.001           # Dimensional total enzyme mass 
𝓒     = 1.0
𝓢     = 1000.0
D_C   = 0.001  # Monomer/polymer diffusivity
D_S   = 0.001  # Substrate diffusivity
Tᵣstar= 0.01  # Release time
ϕ     = 0.5

dims  = fill(Ngrid, nSpatialDims+1)
derivedParams = derivedParameters(Ω, Ωperp, N, k_Cd, k_Ca, k_Sd, k_Sa, k₁, k₂, k₃, k₄, 𝓒, 𝓢, 𝓔, D_C, D_S, Tᵣstar; checks=true)
@unpack L₀, E₀, h₀, C_b, S_b, δ_C, δ_S, α_C, α_S, C₀, S₀, Tᵣ, K₂, K₃, K₄, σ, ϵ, 𝓓, β = derivedParams

#%%

nSpatialDims = 1
Tᵣ = 30.0
K₂ = 1.0
K₄ = 0.0001
α_C = 1.0
𝓓 = 1.0
β = 0.1
Ngrid = 51
dims  = fill(Ngrid, nSpatialDims+1)

sol = glycosylationAnyD(dims, K₂, K₄, Tᵣ, α_C, 𝓓, β, thickness="GRF")

println("finished sim")

#%%

# Create directory for run data labelled with current time.
paramsName = @savename nSpatialDims K₂ K₃ K₄ α_C δ_C σ N β 𝓓 Tᵣ h₀ Ωperp 𝓒
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
    uSlices = [reshape(u, dims...)[:,:,dims[3]÷2] for u in sol.u]
    uSlicesReshaped = [reshape(u, prod(dims[Not(3)])) for u in uSlices]
    concentrationSurfaceMovie(uSlicesReshaped, sol.t, dims[1:2]; subFolder=subFolder, folderName=folderName)
end
