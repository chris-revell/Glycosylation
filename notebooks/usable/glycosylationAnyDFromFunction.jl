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

@from "$(srcdir("Glycosylation.jl"))" using Glycosylation
@from "$(srcdir("Visualise.jl"))" using Visualise
@from "$(srcdir("UsefulFunctions.jl"))" using UsefulFunctions

nSpatialDims = 2

Ωperp = 100.0  # Lumen footprint area
N     = 100         # Maximum polymer length 
k_Cd  = 0.9 # Complex desorption rate
k_Ca  = 1.1 # Complex adsorption rate
k_Sd  = 0.9 # Substrate desorption rate
k_Sa  = 1.1 # Substrate adsorption rate
k₁    = 1.0   # Complex formation forward reaction rate 
k₂    = 0.6   # Complex dissociation reverse reaction rate 
k₃    = 1.1   # Product formation
k₄    = 0.6  # Product dissociation 
E_0   = 1.0
𝓒     = 100.0
𝓢     = 100.0
D_C   = 1.0  # Monomer/polymer diffusivity
D_S   = 1.0  # Substrate diffusivity
Tᵣstar= 100.0  # Release time
ϕ     = 0.5

Nghost   = 1           # Number of ghost points on each side of the domain 
Ngrid = 101

# Create directory for run data labelled with current time.
paramsName = @savename cisternaSeriesID k₂ k₃ k₄ nSpatialDims
folderName = "$(paramsName)_$(Dates.format(Dates.now(),"yy-mm-dd-HH-MM-SS"))"
# Create frames subdirectory to store system state at each output time
subFolder = ""
mkpath(datadir("sims",subFolder,folderName))

sol = glycosylationAnyD(nSpatialDims, Ngrid, Nghost, Ωperp, N, k_Cd, k_Ca, k_Sd, k_Sa, k₁, k₂, k₃, k₄, E_0, 𝓒, 𝓢, D_C, D_S, Tᵣstar, ϕ)

println("finished sim")

if nSpatialDims==1
    concentrationSurfaceMovie(sol.u, sol.t, xs, νs, dimsReal, Nghost, ghostVertexMaskVec; subFolder=subFolder, folderName=folderName)
else
    # uMats = [reshape(u, dimsPlus...) for u in sol.u]
    # uSlices = [selectdim(u, 3, dimsPlus[3]÷2) for u in uMats]
    # concentrationSurfaceMovie([reshape(u, Nνplus*Nxplus) for u in uSlices], sol.t, xs, νs, dimsReal, Nghost, ghostVertexMaskVec; subFolder=subFolder, folderName=folderName)
    spaceIntegralOver_ν_Movie(sol.u, sol.t, xs, νs, dimsReal, Nghost, W, ghostVertexMaskVec; subFolder=subFolder, folderName=folderName)
end

