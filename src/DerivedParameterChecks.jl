
module DerivedParameterChecks

# h₀    = 0.1
# Ωperp = 100.0  # Lumen footprint area
# N     = 100         # Maximum polymer length 
# k_Cd  = 200.0 # Complex desorption rate
# k_Ca  = 1.0 # Complex adsorption rate
# k_Sd  = 200.0 # Substrate desorption rate
# k_Sa  = 1.1 # Substrate adsorption rate
# k₁    = 1.0   # Complex formation forward reaction rate 
# k₂    = 0.1   # Complex dissociation reverse reaction rate 
# k₃    = 1.0   # Product formation
# k₄    = 1.0  # Product dissociation 
# E_0   = 0.001
# 𝓒     = 100.0
# 𝓢     = 1000.0
# D_C   = 0.01  # Monomer/polymer diffusivity
# D_S   = 0.01  # Substrate diffusivity
# Tᵣstar= 10.0  # Release time
# ϕ     = 0.5

function derivedParameterChecks(h₀, Ωperp, N, k_Cd, k_Ca, k_Sd, k_Sa, k₁, k₂, k₃, k₄, E_0, 𝓒, 𝓢, D_C, D_S, Tᵣstar, ϕ)

    𝓔    = 2*Ωperp*E_0   # Total enzyme mass
    δ_C  = π*D_C/(k₁*𝓔)
    δ_S  = π*D_S/(k₁*𝓔)
    Tᵣ   = k₁*𝓔*Tᵣstar/(2*Ωperp)
    Ω    = h₀*Ωperp         # Lumen volume
    α_C  = (k_Cd*Ω)/(2*k_Ca*Ωperp) # Balance of complex in bulk to complex on membrane    
    α_S  = (k_Sd*Ω)/(2*k_Sa*Ωperp) # Balance of substrate in bulk to substrate on membrane 
    C_b  = 𝓒/Ω 
    S_b  = 𝓢/Ω 
    C_0  = C_b*h₀/(2*(1+α_C))      # Early surface monomer concentration
    S_0  = S_b*h₀/(2*(1+α_S))      # Early surface substrate concentration 
    K₂   = (k₂/(k₁*C_b))*((2*k_Ca*Ωperp + k_Cd*Ω)/(k_Ca*Ω)) # Non-dimensionalised complex formation net reaction rate
    K₃   = k₃/k₁    # Non-dimensionalised product formation rate
    K₄   = k₄/k₁    # Non-dimensionalised prodict dissociation rate
    σ    = (k_Sa*S_b*(2*k_Ca*Ωperp + k_Cd*Ω)) / (k_Ca*C_b*(2*k_Sa*Ωperp + k_Sd*Ω))
    ϵ    = 𝓔*(2*k_Ca*Ωperp + k_Cd*Ω) / (2*k_Ca*C_b*Ωperp)
    𝓓    = α_C*δ_C*N^2*(K₂ + σ*K₃)
    β    = N*(σ*K₃ - K₂*K₄)
    L₀  = sqrt(π)*Ω / (Ωperp)^(1.5) # Mean radius 


    println("Small aspect ratio: Ω² << Ω⟂³min(1, D_C/k₁𝓔, D_S/k₁𝓔)")
    # println("Ω² = $(Ω^2), Ω⟂³min(1, D_C/k₁𝓔, D_S/k₁𝓔) = $(Ωperp^3*minimum([1.0,D_C/k₁*𝓔,D_S/k₁*𝓔]))")
    printstyled("$(Ω^2 < Ωperp^3*minimum([1.0,D_C/k₁*𝓔,D_S/k₁*𝓔]))"; color = (Ω^2 < Ωperp^3*minimum([1.0,D_C/k₁*𝓔,D_S/k₁*𝓔]) ? :green : :red))
    println("")

    println("Strong exchange kinetics: D_C*Ωperp << k_Ca*Ω, D_S*Ωperp << k_Sa*Ω") 
    # println("D_C*Ωperp = $(D_C*Ωperp), k_Ca*Ω = $(k_Ca*Ω)")
    printstyled("$(D_C*Ωperp<k_Ca*Ω) "; color = (D_C*Ωperp<k_Ca*Ω ? :green : :red))
    printstyled("$(D_S*Ωperp<k_Sa*Ω)"; color = (D_S*Ωperp<k_Sa*Ω ? :green : :red))
    println("")

    println("Limited enzyme: ϵ << 1 ")
    # println("ϵ = $(ϵ) ")
    printstyled("$(ϵ<1)"; color = (ϵ<1 ? :green : :red))
    println("")

    println("Abundant substrate: σ >> 1")
    # println("σ = $(σ)")
    printstyled("$(σ>1)"; color = (σ>1 ? :green : :red))
    println("")

    println("Slow adsorption of cargo: α_C >> 1") 
    # println("α_C=$(α_C)")
    printstyled("$(α_C>1)"; color = (α_C>1 ? :green : :red))
    println("")

    println("CONFIRM WHERE THIS ONE COMES FROM")
    println("Abundant substrate: k₂k₄k_Sd < S_bk₁k₃k_Sa")
    # println("k₂k₄k_Sd = $(k₂*k₄*k_Sd), S_bk₁k₃k_Sa = $(S_b*k₁*k₃*k_Sa)")
    printstyled("$(k₂*k₄*k_Sd < S_b*k₁*k₃*k_Sa)"; color = (k₂*k₄*k_Sd < S_b*k₁*k₃*k_Sa ? :green : :red))
    println("")

    println("Balanced production: k₄ ∼ k₁")
    # println("k₄ = $(k₄) ∼ k₁ = $(k₁) ")
    printstyled("$(isapprox(k₄, k₁, rtol = 0.05))"; color = (isapprox(k₄, k₁, rtol = 0.05) ? :green : :red))
    println("")

    println("Balanced production: k₁*k_Ca*C_b*(2*k_Sa*Ωperp + k_Sd*Ω) ∼ k₃*k_Sa*S_b*(2*k_Ca*Ωperp + k_Cd*Ω) ")
    println("k₁*k_Ca*C_b*(2*k_Sa*Ωperp + k_Sd*Ω) = $(k₁*k_Ca*C_b*(2*k_Sa*Ωperp + k_Sd*Ω)), k₃*k_Sa*S_b*(2*k_Ca*Ωperp + k_Cd*Ω) = $(k₃*k_Sa*S_b*(2*k_Ca*Ωperp + k_Cd*Ω))")
    printstyled("$(isapprox(k₁*k_Ca*C_b*(2*k_Sa*Ωperp + k_Sd*Ω), k₃*k_Sa*S_b*(2*k_Ca*Ωperp + k_Cd*Ω), rtol = 0.05))"; color = (isapprox(k₁*k_Ca*C_b*(2*k_Sa*Ωperp + k_Sd*Ω), k₃*k_Sa*S_b*(2*k_Ca*Ωperp + k_Cd*Ω), rtol = 0.05) ? :green : :red))
    println("")

    println("Adequate adsorbed substrate: 2k₂k₄k_SaΩperp < (S_bk₁k₃k_Sa - k₂k₄k_Sd)Ω") 
    println("2k₂k₄k_SaΩperp = $(2*k₂*k₄*k_Sa*Ωperp), (S_bk₁k₃k_Sa - k₂k₄k_Sd)Ω=$((S_b*k₁*k₃*k_Sa - k₂*k₄*k_Sd)*Ω)")
    printstyled("$(2*k₂*k₄*k_Sa*Ωperp < (S_b*k₁*k₃*k_Sa - k₂*k₄*k_Sd)*Ω)"; color = (2*k₂*k₄*k_Sa*Ωperp < (S_b*k₁*k₃*k_Sa - k₂*k₄*k_Sd)*Ω ? :green : :red))
    println("")


    println("σK₃ ∼ K₄ ∼ 1")
    println("σK3=$(σ*K₃), K4=$(K₄)")
    printstyled("$(isapprox(σ*K₃, 1.0, rtol=0.1)), "; color = (isapprox(σ*K₃, 1.0, rtol=0.1) ? :green : :red))
    printstyled("$(isapprox(K₄, 1.0, rtol=0.1))"; color = (isapprox(K₄, 1.0, rtol=0.1) ? :green : :red))
    println("")


    return Dict("𝓔"=>𝓔, "K₃"=>K₃, "K₄"=>K₄, "δ_C"=>δ_C, "δ_S"=>δ_S, "Tᵣ"=>Tᵣ, "Ω"=>Ω, "α_C"=>α_C, "α_S"=>α_S, "C_b"=>C_b, "S_b"=>S_b, "C_0"=>C_0, "S_0"=>S_0, "K₂"=>K₂, "σ"=>σ, "ϵ"=>ϵ, "𝓓"=>𝓓, "β"=>β, "K₂"=>K₂, "L₀"=>L₀)
end 

export derivedParameterChecks

end