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

@from "$(srcdir("Glycosylation.jl"))" using Glycosylation
@from "$(srcdir("Visualise.jl"))" using Visualise
@from "$(srcdir("UsefulFunctions.jl"))" using UsefulFunctions
@from "$(srcdir("MakeWeightMatrices.jl"))" using MakeWeightMatrices
@from "$(srcdir("DerivedParameters.jl"))" using DerivedParameters


nSpatialDims = 2

h₀ = 0.0001

# Ωperp = 1.0  # Lumen footprint area
# N     = 100         # Maximum polymer length 
# k_Cd  = 3000.0 # Complex desorption rate
# k_Ca  = 0.01 # Complex adsorption rate
# k_Sd  = 1.0 # Substrate desorption rate
# k_Sa  = 1.0 # Substrate adsorption rate
# k₁    = 2.0   # Complex formation forward reaction rate 
# k₂    = 0.01   # Complex dissociation reverse reaction rate 
# k₃    = 0.01   # Product formation
# k₄    = 2.0  # Product dissociation 
# E_0   = 0.01
# 𝓒     = 1.0
# 𝓢     = 1000.0
# D_C   = 0.000001  # Monomer/polymer diffusivity
# D_S   = 0.000001  # Substrate diffusivity
# Tᵣstar= 0.1  # Release time
# ϕ     = 0.5

# Ngrid = 11
nSpatialDims = 1

Ωperp = 1.0  # Lumen footprint area
N     = 100         # Maximum polymer length 
k_Cd  = 3000.0 # Complex desorption rate
k_Ca  = 0.01 # Complex adsorption rate
k_Sd  = 1.0 # Substrate desorption rate
k_Sa  = 1.0 # Substrate adsorption rate
k₁    = 2.0   # Complex formation forward reaction rate 
k₂    = 0.01   # Complex dissociation reverse reaction rate 
k₃    = 0.01   # Product formation
k₄    = 2.0  # Product dissociation 
E_0   = 0.01
𝓒     = 1.0
𝓢     = 1000.0
D_C   = 0.000001  # Monomer/polymer diffusivity
D_S   = 0.000001  # Substrate diffusivity
Tᵣstar= 0.001  # Release time
ϕ     = 0.5

Ngrid = 101
nSpatialDims == 1 ? dims  = [Ngrid, Ngrid] : dims  = [Ngrid, Ngrid, Ngrid]
derivedParams = derivedParameters(h₀, Ωperp, N, k_Cd, k_Ca, k_Sd, k_Sa, k₁, k₂, k₃, k₄, E_0, 𝓒, 𝓢, D_C, D_S, Tᵣstar; checks=true)
@unpack 𝓔, K₃, K₄, δ_C, δ_S, Tᵣ, Ω, α_C, α_S, C_b, S_b, C_0, S_0, K₂, σ, ϵ, 𝓓, β, K₂, L₀ = derivedParams
#%%

xMax = (Ωperp)^(1/nSpatialDims)
xs   = collect(range(0.0, xMax, dims[2]))
dx   = xs[2]-xs[1]
if nSpatialDims > 1 
    yMax = xMax
    ys   = collect(range(0.0, yMax, dims[3]))
    dy   = ys[2]-ys[1]
end
νMax = 1.0
νs   = collect(range(0.0, νMax, dims[1])) # Positions of discretised vertices in polymerisation space 
dν   = νs[2]-νs[1]
nSpatialDims == 1 ? spacing  = [dν, dx] : spacing  = [dν, dx, dy]

# K₂ = 1.0
# K₃ = 2.0
# K₄ = 1.0  
# α_C = 100.0
# δ_C = 0.01
# σ = 10.0
# N = 100
# β = N*(σ*K₃ - K₂*K₄)
# 𝓓 = α_C*δ_C*N^2*(K₂+σ*K₃)
# Tᵣ = 0.002
# h₀ = 0.01
# Ωperp = 1.0
# 𝓒 = 10.0

mat_h = h₀.*ones(fill(Ngrid, nSpatialDims+1)...)

sol = glycosylationAnyD(mat_h, dims, Ωperp, 𝓒, K₂, K₄, Tᵣ, α_C , 𝓓, β)

println("finished sim")

#%%

# Create directory for run data labelled with current time.
paramsName = @savename nSpatialDims K₂ K₃ K₄ α_C δ_C σ N β 𝓓 Tᵣ h₀ Ωperp 𝓒
folderName = "$(Dates.format(Dates.now(),"yy-mm-dd-HH-MM-SS"))_$(paramsName)"
# Create frames subdirectory to store system state at each output time
subFolder = ""
mkpath(datadir("sims",subFolder,folderName))

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

