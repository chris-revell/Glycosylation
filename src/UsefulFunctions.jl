
module UsefulFunctions

using LinearAlgebra
using SparseArrays
using UnPack
using Statistics

# Integral of h*C over space 
# hᵥ here is dimensionless thickness varying around mean of 1.0 and vertex weights are dimensinless
# Eq 2.47
function M_tilde(u, W, dims, dν, hᵥ)
    uInternal = reshape(W*hᵥ*u, dims...)
    return sum(uInternal, dims=(2:length(dims)))./dν
end

function T_r_star(T̃ᵣ, N, 𝓔, Ω, Ωperp, C_b, S_b, k₁, k₂, k₃, k_Ca, k_Cd, k_Sa, k_Sd)
    T_r_star =  ((2.0*Ωperp*N^2*T̃ᵣ)/𝓔) *
        ((2.0*k_Sa*Ωperp + k_Sd*Ω) / (2.0*k_Ca*Ωperp + k_Cd*Ω)) *    
        (k₁*k_Ca*C_b*Ω) / (k₁*k₂*(2.0*k_Sa*Ωperp + k_Sd*Ω) + k₃*k_Sa*S_b*Ω)
    return T_r_star
end

# Dimensional bulk functional mass integrated over space and polymerisation 
function M_star_ϕ(u, W, dims, dν, hᵥ, α_C, C_b, Ω, ϕ)
    M̃ = M_tilde(u, W, dims, dν, hᵥ)
    Mϕ = dν*sum(M̃[ceil(Int, ϕ*dims[1]) : dims[1]])
    prefactor = α_C*C_b*Ω/(π*(1+α_C))
    return prefactor*Mϕ
end

function P_star(u, W, dims, dν, hᵥ, α_C, C_b, Ω, ϕ, Ωperp, k₁, 𝓔, Tᵣ)
    # Tstar = (2.0*Ωperp*N^2*Tᵣ / 𝓔) * ((2.0*k_Sa*Ωperp + k_Sd*Ω)/(2.0*k_Ca*Ωperp + k_Cd*Ω)) * ((k₁*k_Ca*C_b*Ω)/ (k₁*k₂*(2.0*k_Sa*Ωperp +k_Sd*Ω) + k₃*k_Sa*S_b*Ω))
    # return M_star_ϕ(u, W, dims, dν, hᵥ, α_C, C_b, Ω, ϕ)/Tstar #Dimensional or non-dimensionalised time?
    return (k₁*𝓔/(2*Ωperp))*M_star_ϕ(u, W, dims, dν, hᵥ, α_C, C_b, Ω, ϕ)/Tᵣ 
end

function 𝓟starUniform(𝓒, 𝓔, 𝓢, ϕ, N, k₁, k₂, K₃, K₄, Ωperp, h₀, h_C, h_S)
    Δ = k₁*𝓒/(2.0*k₂*Ωperp)
    return (𝓒/(4*ϕ))*(k₁*𝓔/(Ωperp*N))* (h₀/((h₀+h_C)*(h₀+h_C*(1+Δ)))) * ((𝓢*Δ*K₃/𝓒 - K₄ - K₄*h₀/h_S)/(1 + 𝓢*Δ*K₃/𝓒 + h₀/h_S) )
end

function homogeneousWidthC(ν̃, K₂, K₄, α_C, β, t)
    Ẽ = K₂/(1+K₂)
    M = 1.0
    ξ = ν̃ - Ẽ*β*t/(1+α_C)
    D = Ẽ*K₂*K₄/(1+α_C)
    return (M/sqrt(4.0*π*D*t))*exp(-ξ^2/(4.0*D*t))
end

export T_r_star
export M_tilde
export M_star_ϕ
export P_star
export 𝓟starUniform
export homogeneousWidthC

end