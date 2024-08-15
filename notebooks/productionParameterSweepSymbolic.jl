using Symbolics
using CairoMakie

# @variables Ω Ωperp N k_Cd k_Ca k_Sd k_Sa k₁ k₂ k₃ k₄ C_b S_b S_0 E_0 D_C D_S Tᵣ⁰ α_C α_S K₂  K₃  K₄  h₀  L₀  C_0 𝓒   σ   𝓢   𝓔   ϵ   β δ_C δ_S 𝓓   Tᵣ  ϕ
@variables ϕ Ω Ωperp N k_Cd k_Ca k_Sd k_Sa k₁ k₂ k₃ k₄ C_b S_b S_0 E_0 D_C D_S Tᵣ⁰ 

α_C = (k_Cd*Ω)/(2*k_Ca*Ωperp) # Balance of complex in bulk to complex on membrane
α_S = (k_Sd*Ω)/(2*k_Sa*Ωperp) # Balance of substrate in bulk to substrate on membrane
K₂  = (k₂/(k₁*C_b))*((2*k_Ca*Ωperp + k_Cd*Ω)/(k_Ca*Ω)) # Non-dimensionalised complex formation net reaction rate
K₃  = k₃/k₁    # Non-dimensionalised product formation rate
K₄  = k₄/k₁    # Non-dimensionalised prodict dissociation rate
h₀  = Ω/Ωperp                   # Mean thickness 
L₀  = sqrt(π)*Ω / (Ωperp)^(1.5) # Mean radius 
C_0 = C_b*h₀/(2*(1+α_C))        # Early surface monomer concentration
𝓒   = C_b*Ω                     # Initial monomer mass
σ   = (k_Sa*S_b*(2*k_Ca*Ωperp + k_Cd*Ω)) / (k_Ca*C_b*(2*k_Sa*Ωperp + k_Sd*Ω))
𝓢   = S_b*Ω                     # Initial substrate mass
𝓔   = 2*E_0*Ωperp               # Total enzyme mass
ϵ   = 𝓔*(2*k_Ca*Ωperp + k_Cd*Ω) / (2*k_Ca*C_b*Ωperp)
β   = N*(σ*K₃ - K₂*K₄)
δ_C = π*D_C/(k₁*𝓔)
δ_S = π*D_S/(k₁*𝓔)
𝓓   = α_C*δ_C*N^2*(K₂ + σ*K₃)
Tᵣ  = (k₁*𝓔*Tᵣ⁰)/(2*Ωperp)

𝓟 = (π/(2*ϕ)) * (α_C*C_b*Ω/(1+α_C)^2) * (k₁*𝓔/(2*Ωperp)) * (K₂/(1+K₂)) * ((σ*K₃ - K₂*K₄)/(N*(K₂+σ*K₃)))

varsVec = [ϕ, Ω, Ωperp, N, k_Cd, k_Ca, k_Sd, k_Sa, k₁, k₂, k₃, k₄, C_b, S_b, S_0, E_0, D_C, D_S]

𝓟′ = Symbolics.gradient(𝓟, varsVec)

𝓟′smp = Symbolics.simplify.(𝓟′)





#%%

valuesDict = Dict(
    ϕ => 0.5,
    Ω => 1.0,
    Ωperp => 100.0,
    N => 100,
    k_Cd => 100.0,
    k_Ca => 1.0,
    k_Sd => 1.0,
    k_Sa => 100.0,
    k₁ => 1.0,
    k₂ => 1.0,
    k₃ => 1.0,
    k₄ => 1.0,
    C_b => 1.0,
    S_b => 100.0,
    S_0 => 1.0,
    E_0 => 0.001/Ωperp,
    D_C => 0.001,
    D_S => 0.01,
    Tᵣ⁰ => 1.0,
)

tst = substitute(𝓟, valuesDict)

Symbolics.solve_for(𝓟′[3]~0, Ωperp; simplify=true, check=true)
# tst = Symbolics.solve_for.(𝓟′smp, varsVec)
# tst = Symbolics.solve_for(𝓟′smp[1], varsVec[1])