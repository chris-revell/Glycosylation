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
@from "$(srcdir("MakeWeightMatrices.jl"))" using MakeWeightMatrices
@from "$(srcdir("DerivedParameterChecks.jl"))" using DerivedParameterChecks


nSpatialDims = 2

h₀ = 0.02

Ωperp = 100.0  # Lumen footprint area
N     = 100         # Maximum polymer length 
k_Cd  = 2000.0 # Complex desorption rate
k_Ca  = 2.0 # Complex adsorption rate
k_Sd  = 20.0 # Substrate desorption rate
k_Sa  = 1.1 # Substrate adsorption rate
k₁    = 1.0   # Complex formation forward reaction rate 
k₂    = 0.02   # Complex dissociation reverse reaction rate 
k₃    = 0.01   # Product formation
k₄    = 1.0  # Product dissociation 
E_0   = 0.001
𝓒     = 100.0
𝓢     = 10000.0
D_C   = 0.01  # Monomer/polymer diffusivity
D_S   = 0.01  # Substrate diffusivity
Tᵣstar= 10000.0  # Release time
ϕ     = 0.5

Nghost= 1           # Number of ghost points on each side of the domain 
Ngrid = 51

xMax = 100.0
xs   = collect(range(0.0, xMax, Ngrid+2*Nghost)) # Positions of discretised vertices in space
mat_h = h₀.*ones(fill(Ngrid+2*Nghost, nSpatialDims+1)...)

Nνplus   = Ngrid+2*Nghost # Number of discretised points including ghost points 
Nxplus   = Ngrid+2*Nghost # Number of discretised points including ghost points
Nyplus   = Ngrid+2*Nghost # Number of discretised points including ghost points
nSpatialDims == 1 ? dimsPlus = [Nνplus, Nxplus] : dimsPlus = [Nνplus, Nxplus, Nyplus]
nSpatialDims == 1 ? dimsReal = [Ngrid, Ngrid] : dimsReal = [Ngrid, Ngrid, Ngrid]
dx   = xs[2]-xs[1]
if nSpatialDims > 1 
    yMax = xMax
    ys   = collect(range(0.0, yMax, Nyplus))
    dy   = ys[2]-ys[1]
end
νMax = 1.0
νs   = collect(range(0.0, νMax, Nνplus)) # Positions of discretised vertices in polymerisation space 
dν   = νs[2]-νs[1]
nSpatialDims == 1 ? spacing  = [dν, dx] : spacing  = [dν, dx, dy]

derivedParameterChecks(h₀, Ωperp, N, k_Cd, k_Ca, k_Sd, k_Sa, k₁, k₂, k₃, k₄, E_0, 𝓒, 𝓢, D_C, D_S, Tᵣstar)

#%%

sol = glycosylationAnyD(xs, mat_h, nSpatialDims, Ngrid, Nghost, Ωperp, N, k_Cd, k_Ca, k_Sd, k_Sa, k₁, k₂, k₃, k₄, E_0, 𝓒, 𝓢, D_C, D_S, Tᵣstar, ϕ)

println("finished sim")

#%%

# Create directory for run data labelled with current time.
paramsName = @savename nSpatialDims Ωperp k_Cd k_Ca k_Sd k_Sa k₁ k₂ k₃ k₄ E_0 𝓒 𝓢 D_C D_S Tᵣstar ϕ
folderName = "$(Dates.format(Dates.now(),"yy-mm-dd-HH-MM-SS"))_$(paramsName)"
# Create frames subdirectory to store system state at each output time
subFolder = ""
mkpath(datadir("sims",subFolder,folderName))

ghostVertexMaskVec = makeGhostVertexMask(dimsPlus)
W = vertexVolumeWeightsMatrix(dimsPlus, spacing)

if nSpatialDims==1
    concentrationSurfaceMovie(sol.u, sol.t, dimsReal, Nghost, ghostVertexMaskVec; subFolder="", folderName=folderName)
else
    # uMats = [reshape(u, dimsPlus...) for u in sol.u]
    # uSlices = [selectdim(u, 3, dimsPlus[3]÷2) for u in uMats]
    # concentrationSurfaceMovie([reshape(u, Nνplus*Nxplus) for u in uSlices], sol.t, xs, νs, dimsReal, Nghost, ghostVertexMaskVec; subFolder=subFolder, folderName=folderName)
    spaceIntegralOver_ν_Movie(sol.u, sol.t, xs, νs, dimsReal, Nghost, W, ghostVertexMaskVec; subFolder=subFolder, folderName=folderName)
end

