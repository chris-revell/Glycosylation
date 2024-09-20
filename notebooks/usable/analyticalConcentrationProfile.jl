
using FromFile
using DrWatson
using CairoMakie
using UnPack 

# @from "$(srcdir("UsefulFunctions.jl"))" using UsefulFunctions
@from "$(srcdir("DerivedParameterChecks.jl"))" using DerivedParameterChecks

function homogeneousWidthC(ν, t̃, h₀, Ωperp, N, k_Cd, k_Ca, k_Sd, k_Sa, k₁, k₂, k₃, k₄, E_0, 𝓒, 𝓢, D_C, D_S, Tᵣstar)
    params = derivedParameterNoChecks(h₀, Ωperp, N, k_Cd, k_Ca, k_Sd, k_Sa, k₁, k₂, k₃, k₄, E_0, 𝓒, 𝓢, D_C, D_S, Tᵣstar)
    @unpack 𝓔, Tᵣ, Ω, α_C, C_b, S_b, K₂, K₃, K₄, σ, ϵ, β = params 
    
    Etilde = K₂/(1+K₂)
    
    p1 = (1+α_C)/(π*Etilde*K₂*K₄*t̃)
    p2 = ν*(1+α_C)-Etilde*β*t̃
    p3 = 4*Etilde*K₂*K₄*(1+α_C)*t̃
    
    return sqrt(p1)*exp(-p2^2/p3)
    # return p1

end

Ωperp = 100.0  # Lumen footprint area
N     = 100         # Maximum polymer length 
k_Cd  = 200.0 # Complex desorption rate
k_Ca  = 0.01 # Complex adsorption rate
k_Sd  = 200.0 # Substrate desorption rate
k_Sa  = 1.1 # Substrate adsorption rate
k₁    = 1.0   # Complex formation forward reaction rate 
k₂    = 0.1   # Complex dissociation reverse reaction rate 
k₃    = 1.0   # Product formation
k₄    = 1.0  # Product dissociation 
E_0   = 0.001
𝓒     = 100.0
𝓢     = 1000.0
D_C   = 0.0001  # Monomer/polymer diffusivity
D_S   = 0.0001  # Substrate diffusivity
Tᵣstar= 100.0  # Release time
ϕ     = 0.5

Nghost= 1           # Number of ghost points on each side of the domain 
Ngrid = 51

xMax = 100.0
xs   = collect(range(0.0, xMax, Ngrid+2*Nghost)) # Positions of discretised vertices in space

# h₀s = collect(0.1:0.1:3.0)
h₀s = collect(0.001:0.02:0.2001)

νs = collect(0.00:0.001:1.0)
t̃ = 100.0

params = derivedParameterChecks(h₀s[1], Ωperp, N, k_Cd, k_Ca, k_Sd, k_Sa, k₁, k₂, k₃, k₄, E_0, 𝓒, 𝓢, D_C, D_S, Tᵣstar)

Cs = homogeneousWidthC.(νs, t̃, h₀s[1], Ωperp, N, k_Cd, k_Ca, k_Sd, k_Sa, k₁, k₂, k₃, k₄, E_0, 𝓒, 𝓢, D_C, D_S, Tᵣstar)

#%%

fig = Figure(size=(500,500))
ax = Axis(fig[1,1])
lines!(νs, Cs)
display(fig)