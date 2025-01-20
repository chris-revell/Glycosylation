
module UsefulFunctions

using LinearAlgebra
using SparseArrays
using UnPack
using Statistics

# Integral of h*C over space 
# hᵥ here is dimensionless thickness varying around mean of 1.0 and vertex weights are dimensinless
# \ref{eq:tildeM}
function M_tilde(u, W, dims, dν, hᵥ)
    uInternal = reshape(W*hᵥ*u, dims...)
    return sum(uInternal, dims=(2:length(dims)))./dν
end

function T_r_star(T̃ᵣ, N, ℰ, Ω, Ωperp, C_b, S_b, k₁, k₂, k₃, k_Ca, k_Cd, k_Sa, k_Sd)
    T_r_star =  ((2.0*Ωperp*N^2*T̃ᵣ)/ℰ) *
        ((2.0*k_Sa*Ωperp + k_Sd*Ω) / (2.0*k_Ca*Ωperp + k_Cd*Ω)) *    
        (k₁*k_Ca*C_b*Ω) / (k₁*k₂*(2.0*k_Sa*Ωperp + k_Sd*Ω) + k₃*k_Sa*S_b*Ω)
    return T_r_star
end

# Dimensional bulk functional mass integrated over space and polymerisation 
# \ref{eq:Mstar}
function M_star_ϕ(u, W, dims, dν, hᵥ, α_C, 𝒞, Ω, ϕ)
    M̃ = M_tilde(u, W, dims, dν, hᵥ)
    Mϕ = dν*sum(M̃[floor(Int64, ϕ*dims[1]) : dims[1]])
    prefactor = α_C*𝒞/(π*(1+α_C))
    return prefactor*Mϕ
end

function P_star(u, W, dims, dν, hᵥ, α_C, C_b, Ω, ϕ, Ωperp, k₁, ℰ, TᵣStar)
    return M_star_ϕ(u, W, dims, dν, hᵥ, α_C, C_b, Ω, ϕ)/TᵣStar 
end

function Pstar₅₀Analytic(h₀, h_C, h_S, k₁, k₂, k₃, k₄, Ωperp, 𝒮, 𝒞, ℰ, N, ϕ) 
    u = h₀/h_C
    λ = h_C/h_S
    ζ = (2*k₂*Ωperp)/(k₃*𝒮)
    γ = (2*k₂*Ωperp)/(k₁*𝒞)
    Δ = 2*k₂*k₄*Ωperp/(k₁*k₃*𝒮)
    F = (u*(1-Δ*(1+λ*u)))/((1+u)*(1+ζ*(1+λ*u)*(1+u+(1/γ))))
    return (k₁*𝒞*ℰ*N/ϕ)*F
end

function homogeneousWidthC(ν̃, K₂, K₄, α_C, β, t̃)
    Ẽ = K₂/(1+K₂)
    M = 1.0
    ξ = ν̃ - Ẽ*β*t̃/(1+α_C)
    D = Ẽ*K₂*K₄/(1+α_C)
    return (M/sqrt(4.0*π*D*t̃))*exp(-ξ^2/(4.0*D*t̃))
end

export T_r_star
export M_tilde
export M_star_ϕ
export P_star
# export 𝓟starUniform
export homogeneousWidthC
export Pstar₅₀Analytic

end

# function 𝓟starUniform(𝒞, ℰ, 𝒮, ϕ, N, k₁, k₂, K₃, K₄, Ωperp, h₀, h_C, h_S)
#     Δ = k₁*𝒞/(2.0*k₂*Ωperp)
#     return (𝒞/(4*ϕ))*(k₁*ℰ/(Ωperp*N))* (h₀/((h₀+h_C)*(h₀+h_C*(1+Δ)))) * ((𝒮*Δ*K₃/𝒞 - K₄ - K₄*h₀/h_S)/(1 + 𝒮*Δ*K₃/𝒞 + h₀/h_S) )
# end
