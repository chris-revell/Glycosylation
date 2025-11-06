include("paramsRaw.jl")

using FromFile
@from "$(srcdir("DerivedParameters.jl"))" using DerivedParameters
derivedParams = derivedParameters(Ω, 𝒜, N, k_Cd, k_Ca, k_Sd, k_Sa, k₁, k₂, k₃, k₄, 𝒞, 𝒮, ℰ, D_C, D_S, Tᵣstar; checks=true)
@unpack L₀, E₀, C_b, S_b, δ_C, δ_S, α_C, α_S, C₀, S₀, Tᵣ, T̃ᵣ, K₂, K₃, K₄, σ, ϵ, 𝒟, β, h_C, h_S, λ, ζ, γ, Δ, F = derivedParams
@show L₀
@show E₀
@show C_b
@show S_b
@show δ_C
@show δ_S
@show α_C
@show α_S
@show C₀
@show S₀
@show Tᵣ
@show T̃ᵣ
@show K₂
@show K₃
@show K₄
@show σ
@show ϵ
@show 𝒟
@show β
@show h_C
@show h_S
@show λ
@show ζ
@show γ
@show Δ
@show F



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