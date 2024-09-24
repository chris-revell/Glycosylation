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
using CairoMakie

@from "$(srcdir("Glycosylation.jl"))" using Glycosylation
@from "$(srcdir("Visualise.jl"))" using Visualise
@from "$(srcdir("UsefulFunctions.jl"))" using UsefulFunctions
@from "$(srcdir("MakeWeightMatrices.jl"))" using MakeWeightMatrices
@from "$(srcdir("DerivedParameters.jl"))" using DerivedParameters

nSpatialDims = 1

Ωperp = 100.0  # Lumen footprint area
N     = 100         # Maximum polymer length 
k_Cd  = 200.0 # Complex desorption rate
k_Ca  = 2.0 # Complex adsorption rate
k_Sd  = 200.0 # Substrate desorption rate
k_Sa  = 1.1 # Substrate adsorption rate
k₁    = 1.0   # Complex formation forward reaction rate 
k₂    = 0.1   # Complex dissociation reverse reaction rate 
k₃    = 1.0   # Product formation
k₄    = 1.0  # Product dissociation 
E_0   = 0.001
𝓒     = 100.0
𝓢     = 1000.0
D_C   = 0.01  # Monomer/polymer diffusivity
D_S   = 0.01  # Substrate diffusivity
Tᵣstar= 100.0  # Release time
ϕ     = 0.5

Ngrid = 51
Nν   = Ngrid
Nx   = Ngrid
Ny   = Ngrid
nSpatialDims == 1 ? dims = [Nν, Nx] : dims = [Nν, Nx, Ny]

xMax = 100.0

xs   = collect(range(0.0, xMax, Ngrid)) # Positions of discretised vertices in space

# h₀s = collect(0.1:0.1:3.0)
h₀s = collect(0.001:0.02:0.2001)

sols = []
hᵥs = []
α_Cs = []
Ωs =[]
C_bs =[]

for h₀ in h₀s
    @show h₀

    mat_h = h₀.*ones(fill(Ngrid, nSpatialDims+1)...)

    sol = glycosylationAnyD(xs, mat_h, nSpatialDims, Ngrid, Ωperp, N, k_Cd, k_Ca, k_Sd, k_Sa, k₁, k₂, k₃, k₄, E_0, 𝓒, 𝓢, D_C, D_S, Tᵣstar, ϕ)
    push!(sols, sol)
    Ω    = h₀*Ωperp
    push!(Ωs, Ω)
    α_C  = (k_Cd*Ω)/(2*k_Ca*Ωperp)
    push!(α_Cs, α_C)
    hᵥ_vec = reshape(mat_h, (Ngrid)^(nSpatialDims+1))
    hᵥ = spdiagm(hᵥ_vec)
    push!(hᵥs, hᵥ)
    C_b  = 𝓒/Ω 
    push!(C_bs, C_b)
end

spacing = [1.0/(Ngrid-1), xs[2]-xs[1]]

W = vertexVolumeWeightsMatrix(dims, spacing)
fig = Figure(size=(500,500))
ax1 = Axis(fig[1,1])
Pstars = Float64[]
for i=1:length(sols)
    push!(Pstars, P_star(sols[i][end], W, [Ngrid, Ngrid], hᵥs[i], ϕ, α_Cs[i], C_bs[i], Ωs[i], spacing[1], Tᵣstar))
end
lines!(ax1, h₀s, Pstars)
ax1.xlabel = "h₀"
ax1.ylabel = L"𝓟^*"

Ps = 𝓟starUniform.(ϕ, 𝓒, 𝓢, E_0, h₀s, Ωperp, k_Ca, k_Cd, k_Sa, k_Sd, k₁, k₂, k₃, k₄, N, Tᵣstar)

ax2 = Axis(fig[2,1])
ax2.xlabel = "h₀"
ax2.ylabel = L"𝓟^*"
ylims!(ax2, (0.0,maximum(Ps)))
xlims!(ax2, (0.0,maximum(h₀s)))
lines!(ax2, h₀s, Ps)

linkxaxes!(ax1, ax2)

display(fig)

# save("simulationPvsh.png",fig)



