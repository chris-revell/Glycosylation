
module UsefulFunctions

using LinearAlgebra
using SparseArrays
using UnPack
using Statistics

# Integrate over ν to find E field in spatial dimensions.
# When state vector u is reshaped to an array with shape dims, assume ν is the first dimension of this array
# Function is agnostic about the whether dims is of length 2 or 3.
function E!(u, dims, Esparse, matE, matFₑ, K₂, dν)
    # Convert state vector to matrix of concentrations (We're calculating enzyme distribution, but using bulk concentration?)
    # cs = selectdim(reshape(u, dims...), 1, 2:(dims[1]-1))
    uMat = reshape(u, dims...)
    integ = dν.*(0.5.*selectdim(uMat, 1, 1) .+ dropdims(sum(selectdim(uMat, 1, 2:dims[1]-1), dims=1), dims=1) .+ 0.5.*selectdim(uMat, 1, dims[1]))
    for slice in eachslice(matE, dims=1)
        slice .= matFₑ.*(K₂./(K₂ .+ integ))
    end
    Esparse[diagind(Esparse)] .= reshape(matE, prod(dims))
    return nothing
end

# Function to update linear operator with new values for E at each iteration in solving the ODE system
function updateOperator!(L, u, p, t)
    # @unpack Part1, Part2, u0, dims, Esparse, matE, matFₑ, K₂, dν = p
    E!(u, p.dims, p.Esparse, p.matE, p.matFₑ, p.K₂, p.dν)
    L .= p.Esparse*p.Part1 .+ p.Part2
end

# Integral of h*C over space 
function M_tilde(u, W, dims, dν, hᵥ)
    uInternal = reshape(W*hᵥ*u, dims...)
    return sum(uInternal, dims=(2:length(dims)))./dν
end

# Dimensional bulk functional mass integrated over space and polymerisation 
function M_star(u, W, dims, dν, hᵥ, ϕ, α_C, C_b)
    Ω = π*mean(hᵥ) # Non-dimensionalised Ωperp is always π
    uInternal = reshape(W*hᵥ*u, dims...)
    M̃ = M_tilde(u, W, dims, dν, hᵥ)
    Mϕ = dν*sum(M̃[round(Int, ϕ*dims[1]) : dims[1]])
    prefactor = α_C*C_b*Ω/(π*(1+α_C))
    return prefactor*Mϕ
end

P_star(u, W, dims, dν, hᵥ, ϕ, α_C, C_b, Tᵣ) = M_star(u, W, dims, dν, hᵥ, ϕ, α_C, C_b)/Tᵣ #Dimensional or non-dimensionalised time?
# P_star(u, W, dims, hᵥ, ϕ, α_C, C_b, Ω, dν, Tᵣstar) = M_star(u, W, dims, dν, hᵥ, ϕ, α_C, C_b)/Tᵣstar #Dimensional or non-dimensionalised time?

function 𝓟starUniform(N, h₀, 𝓒, ϕ, E_0, C_b, S_b, Tᵣ, α_C, K₂, K₃, K₄, σ)
# function 𝓟starUniform(ϕ, 𝓒, 𝓢, E_0, h₀, Ωperp, k_Ca, k_Cd, k_Sa, k_Sd, k₁, k₂, k₃, k₄, N, Tᵣstar)
#     𝓔    = 2*Ωperp*E_0
#     Ω    = h₀*Ωperp
#     C_b  = 𝓒/Ω
#     S_b  = 𝓢/Ω
#     Tᵣ   = k₁*𝓔*Tᵣstar/(2*Ωperp)
#     α_C  = (k_Cd*Ω)/(2*k_Ca*Ωperp)
#     K₂   = (k₂/(k₁*C_b))*((2*k_Ca*Ωperp + k_Cd*Ω)/(k_Ca*Ω))
#     K₃   = k₃/k₁
#     K₄   = k₄/k₁
#     σ    = (k_Sa*S_b*(2*k_Ca*Ωperp + k_Cd*Ω)) / (k_Ca*C_b*(2*k_Sa*Ωperp + k_Sd*Ω))
      k₁ = 1.0
      𝓔    = 2*π*E_0
      Ω    = h₀*π
      Ωperp = π
    return π/(2*ϕ) * (α_C*𝓒)/((1+α_C)^2) * (k₁*𝓔)/(2*Ωperp) * K₂/(1+K₂) * (σ*K₃-K₂*K₄)/(N*(K₂+σ*K₃)) * (1/Tᵣ)
end

function homogeneousWidthC(ν̃, K₂, K₄, α_C, β, t)
    Ẽ = K₂/(1+K₂)
    M = 1.0
    ξ = ν̃ - Ẽ*β*t/(1+α_C)
    D = Ẽ*K₂*K₄/(1+α_C)
    return (M/sqrt(4.0*π*D*t))*exp(-ξ^2/(4.0*D*t))
end

export M_tilde
export M_star
export P_star
export E!
export updateOperator!
export 𝓟starUniform
export homogeneousWidthC

end