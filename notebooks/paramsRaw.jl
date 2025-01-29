# h₀ = 1.0
# Ωperp = 10000    # Dimensional lumen footprint area
# Ω     = h₀*Ωperp      # Dimensional lumen volume 
# N     = 100     # Maximum polymer length 
# k_Cd  = 1.0 # Complex desorption rate
# k_Ca  = 0.1 # Complex adsorption rate
# k_Sd  = 1.0 # Substrate desorption rate
# k_Sa  = 1.0 # Substrate adsorption rate
# k₁    = 1.0   # Complex formation forward reaction rate 
# k₂    = 0.025   # Complex dissociation reverse reaction rate 
# k₃    = 0.025   # Product formation
# k₄    = 1.0  # Product dissociation 
# 𝒞     = 10000.0
# 𝒮     = 100000.0
# ℰ     = 0.0000001
# D_C   = 0.00000000002  # Monomer/polymer diffusivity
# D_S   = 0.00000000002  # Substrate diffusivity
# Tᵣstar= 1000000000000000.0  # Release time
# ϕ     = 0.5


h₀ = 1.0
Ωperp = 10000    # Dimensional lumen footprint area
Ω     = h₀*Ωperp      # Dimensional lumen volume 
N     = 100     # Maximum polymer length 
k_Cd  = 0.6 # Complex desorption rate
k_Ca  = 0.1 # Complex adsorption rate
k_Sd  = 0.6 # Substrate desorption rate
k_Sa  = 0.1 # Substrate adsorption rate
k₁    = 1.0   # Complex formation forward reaction rate 
k₂    = 0.02   # Complex dissociation reverse reaction rate 
k₃    = 0.01   # Product formation
k₄    = 1.0  # Product dissociation 
𝒞     = 10000.0
𝒮     = 1000000.0
ℰ     = 0.0000001
D_C   = 0.000000001  # Monomer/polymer diffusivity
D_S   = 0.000000001  # Substrate diffusivity
Tᵣstar= 5000000000000000.0  # Release time
ϕ     = 0.5

using FromFile
@from "$(srcdir("DerivedParameters.jl"))" using DerivedParameters
println(""); println(""); println("")
derivedParams = derivedParameters(Ω, Ωperp, N, k_Cd, k_Ca, k_Sd, k_Sa, k₁, k₂, k₃, k₄, 𝒞, 𝒮, ℰ, D_C, D_S, Tᵣstar; checks=true)
@show h_C = 2*k_Ca/k_Cd
@show h_S = 2*k_Sa/k_Sd
@unpack L₀, E₀, C_b, S_b, δ_C, δ_S, α_C, α_S, C₀, S₀, Tᵣ, T̃ᵣ, K₂, K₃, K₄, σ, ϵ, 𝒟, β, h_C, h_S = derivedParams

