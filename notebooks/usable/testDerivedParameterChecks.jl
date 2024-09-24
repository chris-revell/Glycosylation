
using FromFile

@from "$(srcdir("DerivedParameters.jl"))" using DerivedParameters

h₀    = 0.1
Ωperp = 100.0  # Lumen footprint area
N     = 100         # Maximum polymer length 
k_Cd  = 200.0 # Complex desorption rate
k_Ca  = 1.0 # Complex adsorption rate
k_Sd  = 200.0 # Substrate desorption rate
k_Sa  = 1.1 # Substrate adsorption rate
k₁    = 1.0   # Complex formation forward reaction rate 
k₂    = 0.1   # Complex dissociation reverse reaction rate 
k₃    = 0.1   # Product formation
k₄    = 1.0  # Product dissociation 
E_0   = 0.001
𝓒     = 100.0
𝓢     = 1000.0
D_C   = 0.01  # Monomer/polymer diffusivity
D_S   = 0.01  # Substrate diffusivity
Tᵣstar= 10.0  # Release time
ϕ     = 0.5

params = derivedParameters(h₀, Ωperp, N, k_Cd, k_Ca, k_Sd, k_Sa, k₁, k₂, k₃, k₄, E_0, 𝓒, 𝓢, D_C, D_S, Tᵣstar; checks=true)
