using CairoMakie
using DrWatson

# Basic parameters: geometry
h₀ = 3.0     # Cisternal thickness 
# Ω = 1.0         # Lumen volume
Ωperp = 100.0  # Lumen footprint area
N = 100         # Maximum polymer length 
# Basic parameters: rate constants
k_Cd = 1000.0 # Complex desorption rate
k_Ca = 1.0 # Complex adsorption rate
k_Sd = 1.0 # Substrate desorption rate
k_Sa = 1000.0 # Substrate adsorption rate
k₁ = 1.0   # Complex formation forward reaction rate 
k₂ = 1.0   # Complex dissociation reverse reaction rate 
k₃ = 1.0   # Product formation
k₄ = 0.01   # Product dissociation 
# Basic parameters: concentrations 
C_b = 1.0  # Initial bulk monomer concentration
S_b = 100.0  # Initial bulk substrate concentration
S_0 = 1.0  # Early surface substrate concentration 
E_0 = 0.001 # Mean enzyme concentration
# Basic parameters: diffusivities
D_C = 0.001  # Monomer/polymer diffusivity
D_S = 0.01  # Substrate diffusivity
# Basic parameters: Timescale 
Tᵣ⁰ = 1.0  # Release time


# Derived quantities: rates
α_C = (k_Cd*Ω)/(2*k_Ca*Ωperp) # Balance of complex in bulk to complex on membrane
α_S = (k_Sd*Ω)/(2*k_Sa*Ωperp) # Balance of substrate in bulk to substrate on membrane
K₂  = (k₂/(k₁*C_b))*((2*k_Ca*Ωperp + k_Cd*Ω)/(k_Ca*Ω)) # Non-dimensionalised complex formation net reaction rate
K₃  = k₃/k₁    # Non-dimensionalised product formation rate
K₄  = k₄/k₁    # Non-dimensionalised prodict dissociation rate
# Derived quantities: masses and concentrations 
# h₀  = Ω/Ωperp                   # Mean thickness 
Ω = Ωperp*h₀
L₀  = sqrt(π)*Ω / (Ωperp)^(1.5) # Mean radius 
C_0 = C_b*h₀/(2*(1+α_C))        # Early surface monomer concentration
𝓒   = C_b*Ω                     # Initial monomer mass
σ   = (k_Sa*S_b*(2*k_Ca*Ωperp + k_Cd*Ω)) / (k_Ca*C_b*(2*k_Sa*Ωperp + k_Sd*Ω))
𝓢   = S_b*Ω                     # Initial substrate mass
𝓔   = 2*E_0*Ωperp               # Total enzyme mass
ϵ   = 𝓔*(2*k_Ca*Ωperp + k_Cd*Ω) / (2*k_Ca*C_b*Ωperp)
β = N*(σ*K₃ - K₂*K₄)
# Derived quantities: diffusivities
δ_C = π*D_C/(k₁*𝓔)
δ_S = π*D_S/(k₁*𝓔)
𝓓   = α_C*δ_C*N^2*(K₂ + σ*K₃)
# Derived quantities: non-dimensionalised time
Tᵣ  = k₁*𝓔*Tᵣ⁰/(2*Ωperp)

ϕ = 0.5


function 𝓟(ϕ, E_0, Ωperp, Ω, k_Sa, k_Sd, S_b, k_Ca, k_Cd, C_b, k₁, k₂, k₃, k₄, N)

    α_C = (k_Cd*Ω)/(2*k_Ca*Ωperp)
    K₂  = (k₂/(k₁*C_b))*((2*k_Ca*Ωperp + k_Cd*Ω)/(k_Ca*Ω))
    K₃  = k₃/k₁
    K₄  = k₄/k₁
    𝓔 = 2*E_0*Ωperp
    σ = (k_Sa*S_b*(2*k_Ca*Ωperp + k_Cd*Ω)) / (k_Ca*C_b*(2*k_Sa*Ωperp + k_Sd*Ω))
    
    return (π/(2*ϕ)) * (α_C*C_b*Ω/(1+α_C)^2) * (k₁*𝓔/(2*Ωperp)) * (K₂/(1+K₂)) * ((σ*K₃ - K₂*K₄)/(N*(K₂+σ*K₃)))
end 


#%%

@show ϕ
@show α_C
@show C_b
@show Ω
@show Ωperp
@show k₁
@show K₂
@show K₃
@show K₄
@show σ
@show N
@show 𝓔


println("Small aspect ratio")
println("Ω² << Ω⟂³min(1, D_C/k₁𝓔, D_S/k₁𝓔)")
println("Ω² = $(Ω^2), Ω⟂³min(1, D_C/k₁𝓔, D_S/k₁𝓔) = $(Ωperp^3*minimum([1.0,D_C/k₁*𝓔,D_S/k₁*𝓔]))")
# println("$(Ω^2 < Ωperp^3*minimum([1.0,D_C/k₁*𝓔,D_S/k₁*𝓔]))")
printstyled("$(Ω^2 < Ωperp^3*minimum([1.0,D_C/k₁*𝓔,D_S/k₁*𝓔]))"; color = (Ω^2 < Ωperp^3*minimum([1.0,D_C/k₁*𝓔,D_S/k₁*𝓔]) ? :green : :red))
println("")

println("Limited enzyme")
println("ϵ << 1 ")
println("ϵ = $(ϵ) ")
# println("$(ϵ<1)")
printstyled("$(ϵ<1)"; color = (ϵ<1 ? :green : :red))
println("")

println("Abundant substrate")
println("σ >> 1")
println("σ = $(σ)")
# println("$(σ>1)")
printstyled("$(σ>1)"; color = (σ>1 ? :green : :red))
println("")

println("Abundant substrate")
println("k₂k₄k_Sd < S_bk₁k₃k_Sa")
println("k₂k₄k_Sd = $(k₂*k₄*k_Sd), S_bk₁k₃k_Sa = $(S_b*k₁*k₃*k_Sa)")
# println("$(k₂*k₄*k_Sd < S_b*k₁*k₃*k_Sa)")
printstyled("$(k₂*k₄*k_Sd < S_b*k₁*k₃*k_Sa)"; color = (k₂*k₄*k_Sd < S_b*k₁*k₃*k_Sa ? :green : :red))
println("")

# println("Balanced production")
# println("k₄ ∼ k₁")
# println("k₄ = $(k₄) ∼ k₁ = $(k₁) ")
# # println("$(isapprox(k₄, k₁, rtol = 0.05))")
# printstyled("$(isapprox(k₄, k₁, rtol = 0.05))"; color = (isapprox(k₄, k₁, rtol = 0.05) ? :green : :red))
# println("")

# println("Balanced production")
# println("k₁*k_Ca*C_b*(2*k_Sa*Ωperp + k_Sd*Ω) ∼ k₃*k_Sa*S_b*(2*k_Ca*Ωperp + k_Cd*Ω) ")
# println("k₁*k_Ca*C_b*(2*k_Sa*Ωperp + k_Sd*Ω) = $(k₁*k_Ca*C_b*(2*k_Sa*Ωperp + k_Sd*Ω)), k₃*k_Sa*S_b*(2*k_Ca*Ωperp + k_Cd*Ω) = $(k₃*k_Sa*S_b*(2*k_Ca*Ωperp + k_Cd*Ω))")
# # println("$(isapprox(k₁*k_Ca*C_b*(2*k_Sa*Ωperp + k_Sd*Ω), k₃*k_Sa*S_b*(2*k_Ca*Ωperp + k_Cd*Ω), rtol = 0.05))")
# printstyled("$(isapprox(k₁*k_Ca*C_b*(2*k_Sa*Ωperp + k_Sd*Ω), k₃*k_Sa*S_b*(2*k_Ca*Ωperp + k_Cd*Ω), rtol = 0.05))"; color = (isapprox(k₁*k_Ca*C_b*(2*k_Sa*Ωperp + k_Sd*Ω), k₃*k_Sa*S_b*(2*k_Ca*Ωperp + k_Cd*Ω), rtol = 0.05) ? :green : :red))
# println("")

println("Strong exchange kinetics")
println("D_C*Ωperp << k_Ca*Ω") 
println("D_C*Ωperp = $(D_C*Ωperp), k_Ca*Ω = $(k_Ca*Ω)")
# println("$(D_C*Ωperp<k_Ca*Ω)")
printstyled("$(D_C*Ωperp<k_Ca*Ω)"; color = (D_C*Ωperp<k_Ca*Ω ? :green : :red))
println("")

println("Strong exchange kinetics")
println("D_S*Ωperp << k_Sa*Ω") 
println("D_S*Ωperp = $(D_S*Ωperp), k_Sa*Ω = $(k_Sa*Ω)")
# println("$(D_S*Ωperp<k_Sa*Ω)")
printstyled("$(D_S*Ωperp<k_Sa*Ω)"; color = (D_S*Ωperp<k_Sa*Ω ? :green : :red))
println("")

println("Adequate adsorbed substrate")
println("2k₂k₄k_SaΩperp < (S_bk₁k₃k_Sa - k₂k₄k_Sd)Ω") 
println("2k₂k₄k_SaΩperp = $(2*k₂*k₄*k_Sa*Ωperp), (S_bk₁k₃k_Sa - k₂k₄k_Sd)Ω=$((S_b*k₁*k₃*k_Sa - k₂*k₄*k_Sd)*Ω)")
# println("$(2*k₂*k₄*k_Sa*Ωperp < (S_b*k₁*k₃*k_Sa - k₂*k₄*k_Sd)*Ω)")
printstyled("$(2*k₂*k₄*k_Sa*Ωperp < (S_b*k₁*k₃*k_Sa - k₂*k₄*k_Sd)*Ω)"; color = (2*k₂*k₄*k_Sa*Ωperp < (S_b*k₁*k₃*k_Sa - k₂*k₄*k_Sd)*Ω ? :green : :red))
println("")


println("Slow adsorption of cargo")
println("α_C >> 1") 
println("α_C=$(α_C)")
# println("$(α_C>1)")
printstyled("$(α_C>1)"; color = (α_C>1 ? :green : :red))
println("")


@show 𝓟(ϕ, E_0, Ωperp, Ω, k_Sa, k_Sd, S_b, k_Ca, k_Cd, C_b, k₁, k₂, k₃, k₄, N)


#%%

# Ωperps = collect(10:10:500)
# Ωperps = collect(10:10:500)

Ωperp = 100.0
Ωs = collect(10:10:500)

# 𝓟s = 𝓟.(ϕ, E_0, Ωperps, Ω, k_Sa, k_Sd, S_b, k_Ca, k_Cd, C_b, k₁, k₂, k₃, k₄, N)
𝓟s = 𝓟.(ϕ, E_0, Ωperp, Ωs, k_Sa, k_Sd, S_b, k_Ca, k_Cd, C_b, k₁, k₂, k₃, k₄, N)


fig = Figure(size=(750,750), fontsize=16)
ax = Axis(fig[1,1], xlabelsize=32, ylabelsize=32)
ax.xlabel = "Ω"
ax.ylabel = "𝓟"
lines!(Ωs, 𝓟s)
display(fig)
paramsName = @savename ϕ E_0 Ωperp k_Sa k_Sd S_b k_Ca k_Cd C_b k₁ k₂ k₃ k₄ N
save(datadir("$(paramsName).png"), fig)
