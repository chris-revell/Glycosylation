 using Symbolics

 @variables α_C α_S k_Ca k_Cd k_Sa k_Sd Ω Ωperp h₀ L₀ k₁ k₂ k₃ k₄ K₂ K₃ K₄ 
 @variables C_b 𝓒 C_0 S_b 𝓢 S_0 σ E_0 𝓔 ϵ β N δ_C D_C δ_S D_S 𝓓 Tᵣ Tᵣstar Δ ϕ 𝓟star
 @variables Cᵥᵥ Cᵥ Cₜ E


 exprs [
    𝓔    ~ 2*Ωperp*E_0
    Tᵣ   ~ k₁*𝓔*Tᵣstar/(2*Ωperp)
    Ω    ~ h₀*Ωperp
    α_C  ~ (k_Cd*Ω)/(2*k_Ca*Ωperp)
    C_b  ~ 𝓒/Ω 
    S_b  ~ 𝓢/Ω 
    K₂   ~ (k₂/(k₁*C_b))*((2*k_Ca*Ωperp + k_Cd*Ω)/(k_Ca*Ω))
    K₃   ~ k₃/k₁
    K₄   ~ k₄/k₁
    σ    ~ (k_Sa*S_b*(2*k_Ca*Ωperp + k_Cd*Ω)) / (k_Ca*C_b*(2*k_Sa*Ωperp + k_Sd*Ω))
    ϵ    ~ 𝓔*(2*k_Ca*Ωperp + k_Cd*Ω) / (2*k_Ca*C_b*Ω*Ωperp)
    β    ~ N*(σ*K₃ - K₂*K₄)
    E    ~ K₂/(π*(1+K₂))
    p1   ~ (1+α_C)/(4*π*E*K₂*K₄*Tᵣ)
    p2   ~ ν̃*(1+α_C)-E*β*Tᵣ
    p3   ~ 4*E*K₂*K₄*(1+α_C)*Tᵣ
 ]

 Cₜ ~ (E*K₂*K₄*Cᵥᵥ - E*β*Cᵥ)/(1+α_C)

 D ~ E*K₂*K₄/(1+α_C)

 u ~ E*β/(1+α_C)

 M ~ 4*π*D*t/sqrt()
    
 
 
 \tilde{E}(\tilde{t})&={\color{red}\frac{1}{\pi}}\left[1+\frac{1}{K_2}\int_0^\infty \tilde{C}(\tilde{t};\nu)\,\mathrm{d}\nu\right]^{-1}