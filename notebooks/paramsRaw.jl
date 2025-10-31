# # In paper
h₀    = 1.0
𝒜     = 10000    # Dimensional lumen footprint area
Ω     = h₀*𝒜      # Dimensional lumen volume 
N     = 100     # Maximum polymer length 
k_Sa  = 1.0 # Substrate adsorption rate
k_Sd  = 1.0 # Substrate desorption rate
k_Ca  = 0.1 # Complex adsorption rate
k_Cd  = 1.0 # Complex desorption rate
k₁    = 1.0   # Complex formation forward reaction rate 
k₂    = 0.025   # Complex dissociation reverse reaction rate 
k₃    = 0.025   # Product formation
k₄    = 1.0  # Product dissociation 
𝒞     = 1e4 #10000.0
𝒮     = 1e5 #100000.0
ℰ     = 1.0
D_C   = 1e-3 #0.01   # Monomer/polymer diffusivity
D_S   = 1e-3 #0.01   # Substrate diffusivity
Tᵣstar= 1e9 #100000000.0  # Release time
ϕ     = 0.5
differencing = "centre"
solver = SSPRK432()
nOutputs = 1000

# For figure 2:
# h₀    = 1.0
# 𝒜     = 10000    # Dimensional lumen footprint area
# Ω     = h₀*𝒜      # Dimensional lumen volume 
# N     = 100     # Maximum polymer length 
# k_Sa  = 0.1     # Substrate adsorption rate
# k_Sd  = 0.6     # Substrate desorption rate
# k_Ca  = 0.1     # Complex adsorption rate
# k_Cd  = 0.6     # Complex desorption rate
# k₁    = 1.0    # Complex formation forward reaction rate 
# k₂    = 0.02   # Complex dissociation reverse reaction rate 
# k₃    = 0.01   # Product formation
# k₄    = 1.0    # Product dissociation 
# 𝒞     = 10000.0
# 𝒮     = 1000000.0
# ℰ     = 0.0000001
# D_C   = 0.000000001  # Monomer/polymer diffusivity
# D_S   = 0.000000001  # Substrate diffusivity
# Tᵣstar= 2000000000000000.0  # Release time
# ϕ     = 0.5

# using FromFile
# @from "$(srcdir("DerivedParameters.jl"))" using DerivedParameters
# derivedParams = derivedParameters(Ω, 𝒜, N, k_Cd, k_Ca, k_Sd, k_Sa, k₁, k₂, k₃, k₄, 𝒞, 𝒮, ℰ, D_C, D_S, Tᵣstar; checks=true)
# @unpack L₀, E₀, C_b, S_b, δ_C, δ_S, α_C, α_S, C₀, S₀, Tᵣ, T̃ᵣ, K₂, K₃, K₄, σ, ϵ, 𝒟, β, h_C, h_S, λ, ζ, γ, Δ, F = derivedParams
# @show α_C

# @show T₀ = 1/(k₁*E₀)
# @show T̃₀ = T₀/((N^2)*(K₂ + σ*K₃))
# @show sqrt(D_C*T̃₀)
# @show h₀/sqrt(D_C*T̃₀)
# @show h₀/sqrt(D_S*T̃₀)
# @show D_S/(h₀*k_Sa)
# @show D_C/(h₀*k_Ca)
# @show S_b*k₁*k₃*k_Sa/(k₂*k₄*k_Sd)
# hcutoff = (2.0*k_Sa/k_Sd)*((𝒮*k₁*k₃)/(2.0*𝒜*k₂*k₄) - 1.0)
# @show hcutoff