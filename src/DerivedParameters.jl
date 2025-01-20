
module DerivedParameters

function derivedParameters(Ω, Ωperp, N, k_Cd, k_Ca, k_Sd, k_Sa, k₁, k₂, k₃, k₄, 𝒞, 𝒮, ℰ, D_C, D_S, Tᵣstar; checks=true)

    L₀   = sqrt(Ωperp/π)       # Dimensional mean cyclindrical radius of cisterna 
    # ℰ    = 2*Ωperp*E₀        # Dimensional total enzyme mass
    E₀  = ℰ/2*Ωperp           # Dimensional mean enzyme concentration
    # Ω    = h₀*Ωperp           # Dimensional lumen volume
    h₀   = Ω/Ωperp             # Dimensional mean lumen thickness
    C_b  = 𝒞/Ω                 # Dimensional initial bulk monomeric cargo concentration
    S_b  = 𝒮/Ω                 # Dimensional initial bulk substrate concentration
    
    δ_C  = π*D_C/(k₁*ℰ)  # Dimensionless diffusivity
    δ_S  = π*D_S/(k₁*ℰ)  # Dimensionless diffusivity

    α_C  = (k_Cd*Ω)/(2*k_Ca*Ωperp) # Dimensionless complex capacitance
    α_S  = (k_Sd*Ω)/(2*k_Sa*Ωperp) # Dimensionless substrate capacitance    

    C₀   = 𝒞/(2*Ωperp*(1+α_C))  # Dimensional Early surface monomer concentration
    S₀   = 𝒮/(2*Ωperp*(1+α_S))  # Dimensional Early surface substrate concentration 
    
    Tᵣ   = k₁*ℰ*Tᵣstar/(2*Ωperp)   # Dimensionless release time 
    
    K₂   = (k₂/(k₁*C_b))*((2*k_Ca*Ωperp + k_Cd*Ω)/(k_Ca*Ω)) # Dimensionless complex formation net reaction rate
    K₃   = k₃/k₁                                            # Dimensionless product formation rate
    K₄   = k₄/k₁                                            # Dimensionless prodict dissociation rate
    # σ    = (k_Sa*S_b*(2*k_Ca*Ωperp + k_Cd*Ω)) / (k_Ca*C_b*(2*k_Sa*Ωperp + k_Sd*Ω))
    σ    = S₀/C₀                                            # Dimensionless substrate/cargo concentration on surface
    # σ    = 𝒮*(1+α_C)/(𝒞*(1+α_S))
    # ϵ    = ℰ*(2*k_Ca*Ωperp + k_Cd*Ω) / (2*k_Ca*C_b*Ω*Ωperp)
    ϵ    = E₀/C₀                                            # Dimensionless enzyme/cargo concentration on surface 
    # ϵ    = ℰ*(1+α_C)/𝒞
    𝒟    = α_C*δ_C*N^2*(K₂ + σ*K₃)    # Dimensionless parameter on diffusion term, derived from combination of other terms
    β    = N*(σ*K₃ - K₂*K₄)           # Dimensionless parameter on advection term, derived from combination of other terms 

    T̃ᵣ   = Tᵣ/((N^2)*(K₂+σ*K₃))

    h_C = 2*k_Ca/k_Cd
    h_S = 2*k_Sa/k_Sd


    # λ = (𝒮/(2*Ωperp))*(k₁*k₃/(k₂*k₄))

    if checks 
        println("Small aspect ratio: Ω² << Ω⟂³min(1, D_C/k₁ℰ, D_S/k₁ℰ)")
        println("Ω² = $(Ω^2), Ω⟂³min(1, D_C/k₁ℰ, D_S/k₁ℰ) = $((Ωperp^3)*minimum([1.0, D_C/(k₁*ℰ), D_S/(k₁*ℰ)]))")
        printstyled("$(Ω^2 < (Ωperp^3)*minimum([1.0, D_C/(k₁*ℰ), D_S/(k₁*ℰ)]))"; 
            color = (Ω^2 < (Ωperp^3)*minimum([1.0, D_C/(k₁*ℰ), D_S/(k₁*ℰ)]) ? :green : :red))
        println("")

        println("Strong exchange kinetics: D_C*Ωperp << k_Ca*Ω, D_S*Ωperp << k_Sa*Ω") 
        println("D_C*Ωperp = $(D_C*Ωperp), k_Ca*Ω = $(k_Ca*Ω)")
        println("D_S*Ωperp = $(D_S*Ωperp), k_Sa*Ω = $(k_Sa*Ω)")
        printstyled("$(D_C*Ωperp<k_Ca*Ω) "; color = (D_C*Ωperp<k_Ca*Ω ? :green : :red))
        printstyled("$(D_S*Ωperp<k_Sa*Ω)"; color = (D_S*Ωperp<k_Sa*Ω ? :green : :red))
        println("")

        println("Limited enzyme: ϵ << 1 ")
        println("ϵ = $(ϵ) ")
        printstyled("$(ϵ<1)"; color = (ϵ<1 ? :green : :red))
        println("")

        println("Abundant substrate: σ >> 1")
        println("σ = $(σ)")
        printstyled("$(σ>1)"; color = (σ>1 ? :green : :red))
        println("")

        println("Slow adsorption of cargo: α_C >> 1") 
        println("α_C=$(α_C)")
        printstyled("$(α_C>1)"; color = (α_C>1 ? :green : :red))
        println("")

        println("CONFIRM WHERE THIS ONE COMES FROM")
        println("Abundant substrate: k₂k₄k_Sd < S_bk₁k₃k_Sa")
        println("k₂k₄k_Sd = $(k₂*k₄*k_Sd), S_bk₁k₃k_Sa = $(S_b*k₁*k₃*k_Sa)")
        printstyled("$(k₂*k₄*k_Sd < S_b*k₁*k₃*k_Sa)"; color = (k₂*k₄*k_Sd < S_b*k₁*k₃*k_Sa ? :green : :red))
        println("")

        println("Balanced production: k₄ ∼ k₁")
        println("k₄ = $(k₄) ∼ k₁ = $(k₁) ")
        printstyled("$(isapprox(k₄, k₁, rtol = 0.05))"; color = (isapprox(k₄, k₁, rtol = 0.05) ? :green : :red))
        println("")

        println("Balanced production: k₁*k_Ca*C_b*(2*k_Sa*Ωperp + k_Sd*Ω) ∼ k₃*k_Sa*S_b*(2*k_Ca*Ωperp + k_Cd*Ω) ")
        println("k₁*k_Ca*C_b*(2*k_Sa*Ωperp + k_Sd*Ω) = $(k₁*k_Ca*C_b*(2*k_Sa*Ωperp + k_Sd*Ω)), k₃*k_Sa*S_b*(2*k_Ca*Ωperp + k_Cd*Ω) = $(k₃*k_Sa*S_b*(2*k_Ca*Ωperp + k_Cd*Ω))")
        printstyled("$(isapprox(k₁*k_Ca*C_b*(2*k_Sa*Ωperp + k_Sd*Ω), k₃*k_Sa*S_b*(2*k_Ca*Ωperp + k_Cd*Ω), rtol = 0.05))"; color = (isapprox(k₁*k_Ca*C_b*(2*k_Sa*Ωperp + k_Sd*Ω), k₃*k_Sa*S_b*(2*k_Ca*Ωperp + k_Cd*Ω), rtol = 0.1) ? :green : :red))
        println("")

        println("Adequate adsorbed substrate: 2k₂k₄k_SaΩperp < (S_bk₁k₃k_Sa - k₂k₄k_Sd)Ω") 
        println("2k₂k₄k_SaΩperp = $(2*k₂*k₄*k_Sa*Ωperp), (S_bk₁k₃k_Sa - k₂k₄k_Sd)Ω=$((S_b*k₁*k₃*k_Sa - k₂*k₄*k_Sd)*Ω)")
        printstyled("$(2*k₂*k₄*k_Sa*Ωperp < (S_b*k₁*k₃*k_Sa - k₂*k₄*k_Sd)*Ω)"; color = (2*k₂*k₄*k_Sa*Ωperp < (S_b*k₁*k₃*k_Sa - k₂*k₄*k_Sd)*Ω ? :green : :red))
        println("")

        println("σK₃ ∼ K₄ ∼ 1")
        println("σK₃ = $(σ*K₃), K₄ = $(K₄)")
        printstyled("$(isapprox(σ*K₃, 1.0, rtol=0.1)), "; color = (isapprox(σ*K₃, 1.0, rtol=0.1) ? :green : :red))
        printstyled("$(isapprox(K₄, 1.0, rtol=0.1))"; color = (isapprox(K₄, 1.0, rtol=0.1) ? :green : :red))
        println("")

        # println("λ > 1")
        # println("λ = $(λ)")
        # printstyled("$(λ>1)"; color = (λ>1) ? :green : :red)
        # println("")

        # println("h₀ < 2k_Sa(λ-1)/k_Sd")
        # println("h₀ = $(h₀), 2k_Sa(λ-1)/k_Sd = $(2.0*k_Sa*(λ-1)/k_Sd)")
        # printstyled("$(h₀<(2.0*k_Sa*(λ-1)/k_Sd))"; color = (h₀<(2.0*k_Sa*(λ-1)/k_Sd)) ? :green : :red)
        # println("")

    end

    return Dict("L₀"=>L₀, "E₀"=>E₀, "C_b"=>C_b, "S_b"=>S_b, "δ_C"=>δ_C, "δ_S"=>δ_S, "α_C"=>α_C, "α_S"=>α_S, "C₀"=>C₀, "S₀"=>S₀, "Tᵣ"=>Tᵣ, "T̃ᵣ"=>T̃ᵣ, "K₂"=>K₂, "K₃"=>K₃, "K₄"=>K₄, "σ"=>σ, "ϵ"=>ϵ, "𝒟"=>𝒟, "β"=>β, "h_C"=>h_C, "h_S"=>h_S) #, "λ"=>λ)
end 

export derivedParameters

end