
module UsefulFunctions

using LinearAlgebra
using SparseArrays
using UnPack

# Integrate over ν to find E field in spatial dimensions.
# When state vector u is reshaped to an array with shape dimsPlus, assume ν is the first dimension of this array
# Function is agnostic about the whether dimsPlus is of length 2 or 3.
function E!(u, dimsPlus, Esparse, matE, matFₑ, K₂, dν)
    # Convert state vector to matrix of concentrations (We're calculating enzyme distribution, but using bulk concentration?)
    # cs = selectdim(reshape(u, dimsPlus...), 1, 2:(dimsPlus[1]-1))
    uMat = reshape(u, dimsPlus...)
    integ = 0.5.*selectdim(uMat, 1, 1) .+ dropdims(sum(selectdim(uMat, 1, 2:dimsPlus[1]-1), dims=1), dims=1) .+ 0.5.*selectdim(uMat, 1, dimsPlus[1])
    for slice in eachslice(matE, dims=1)
        slice .= matFₑ.*(K₂./(K₂ .+ integ))
    end
    Esparse[diagind(Esparse)] .= reshape(matE, prod(dimsPlus))
    return nothing
end

# Function to update linear operator with new values for E at each iteration in solving the ODE system
function updateOperator!(L, u, p, t)
    @unpack L1, L2, u0, dimsPlus, Esparse, matE, matFₑ, K₂, dν = p
    E!(u, dimsPlus, Esparse, matE, matFₑ, K₂, dν)
    L .= Esparse*L1 .+ L2
end

# Integral of h*C over space 
function M_tilde(u, W, ghostVertexMaskVec, dims, dν, hᵥ)
    uInternal = reshape((W*hᵥ*u)[ghostVertexMaskVec], dims...)
    M_tilde = sum(uInternal, dims=(2:length(dims)))
    return M_tilde./dν
end

# Dimensional bulk functional mass integrated over space and polymerisation 
function M_star(u, W, ghostVertexMaskVec, dims, hᵥ, ϕ, α_C, C_b, Ω, dν)
    uInternal = reshape((W*hᵥ*u)[ghostVertexMaskVec], dims...)
    M̃ = M_tilde(u, W, ghostVertexMaskVec, dims, dν, hᵥ)
    Mϕ = dν*sum(M̃[round(Int, ϕ*dims[1]) : dims[1]])
    prefactor = α_C*C_b*Ω/(π*(1+α_C))
    return prefactor*Mϕ
end

P_star(u, W, ghostVertexMaskVec, dims, hᵥ, ϕ, α_C, C_b, Ω, dν, Tᵣstar) = M_star(u, W, ghostVertexMaskVec, dims, hᵥ, ϕ, α_C, C_b, Ω, dν)/Tᵣstar

function 𝓟starUniform(ϕ, 𝓒, 𝓢, E_0, h₀, Ωperp, k_Ca, k_Cd, k_Sa, k_Sd, k₁, k₂, k₃, k₄, N, Tᵣstar)
    𝓔    = 2*Ωperp*E_0
    Ω    = h₀*Ωperp
    C_b  = 𝓒/Ω
    S_b  = 𝓢/Ω
    Tᵣ   = k₁*𝓔*Tᵣstar/(2*Ωperp)
    α_C  = (k_Cd*Ω)/(2*k_Ca*Ωperp)
    K₂   = (k₂/(k₁*C_b))*((2*k_Ca*Ωperp + k_Cd*Ω)/(k_Ca*Ω))
    K₃   = k₃/k₁
    K₄   = k₄/k₁
    σ    = (k_Sa*S_b*(2*k_Ca*Ωperp + k_Cd*Ω)) / (k_Ca*C_b*(2*k_Sa*Ωperp + k_Sd*Ω))
    return π/(2*ϕ) * (α_C*𝓒)/((1+α_C)^2) * (k₁*𝓔)/(2*Ωperp) * K₂/(1+K₂) * (σ*K₃-K₂*K₄)/(N*(K₂+σ*K₃)) * (1/Tᵣ)
end

function homogeneousWidthC(ν̃, t̃, h₀, 𝓒, k_Ca, k_Cd, k_Sa, k_Sd, k₁, k₂, k₃, k₄, Ωperp, E_0, Tᵣstar)
    𝓔    = 2*Ωperp*E_0
    Tᵣ   = k₁*𝓔*Tᵣstar/(2*Ωperp)
    Ω    = h₀*Ωperp
    α_C  = (k_Cd*Ω)/(2*k_Ca*Ωperp)
    C_b  = 𝓒/Ω 
    S_b  = 𝓢/Ω 
    K₂   = (k₂/(k₁*C_b))*((2*k_Ca*Ωperp + k_Cd*Ω)/(k_Ca*Ω))
    K₃   = k₃/k₁
    K₄   = k₄/k₁
    σ    = (k_Sa*S_b*(2*k_Ca*Ωperp + k_Cd*Ω)) / (k_Ca*C_b*(2*k_Sa*Ωperp + k_Sd*Ω))
    ϵ    = 𝓔*(2*k_Ca*Ωperp + k_Cd*Ω) / (2*k_Ca*C_b*Ω*Ωperp)
    β    = N*(σ*K₃ - K₂*K₄)
    Etilde = K₂/(1+K₂)
    p1 = (1+α_C)/(π*Etilde*K₂*K₄*t̃)
    p2 = ν̃*(1+α_C)-Etilde*β*t̃
    p3 = 4*Etilde*K₂*K₄*(1+α_C)*t̃
    return sqrt(p1)*exp(-p2^2/p3)
end

export M_tilde
export M_star
export P_star
export E!
export updateOperator!
export 𝓟starUniform

end