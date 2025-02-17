
module UsefulFunctions

using LinearAlgebra
using SparseArrays
using UnPack
using Statistics
using SpecialFunctions

# Integral of h*C over space 
# hᵥ here is dimensionless thickness varying around mean of 1.0 and vertex weights are dimensinless
# \ref{eq:tildeM}
function M̃(u, W, dims, dν, hᵥ)
    uInternal = reshape(W*hᵥ*u, dims...)
    return sum(uInternal, dims=(2:length(dims)))./dν
end

function M̃ϕ2(u, W, dims, dν, hᵥ, ϕ; thresh="floor")
    # M̃local = M̃(u, W, dims, dν, hᵥ)
    # M̃ϕ = dν*sum(M̃local[floor(Int64, ϕ*dims[1]) : dims[1]])

    # uInternal = reshape(W*hᵥ*u, dims...)
    # M̃ϕ = sum(selectdim(uInternal, 1, ceil(Int64, ϕ*dims[1]):dims[1]))

    # uInternal = reshape(p1.W*p1.hᵥ*u, p1.dims...)
    # dx = sqrt(π)/(dims[2]-1)
    # tmp = 0.5*dx.*sum(uInternal[:, 1:end-1].+uInternal[:, 2:end], dims=2)
    # tmp = tmp[:,1]
    # dν = 1.0/(dims[1]-1)
    # M̃ϕ = 0.5*dν*sum(tmp[1:end-1].+tmp[2:end])
    
    uInternal = reshape(hᵥ*u, dims...)
    dx = sqrt(π)/(dims[2]-1)
    yMax = sqrt(π)
    if thresh == "floor"
        tmp = 0.5*yMax*dx.*sum(uInternal[floor(Int64, ϕ*dims[1]):dims[1], 1:end-1].+uInternal[floor(Int64, ϕ*dims[1]):dims[1], 2:end], dims=2)
    else 
        tmp = 0.5*yMax*dx.*sum(uInternal[ceil(Int64, ϕ*dims[1]):dims[1], 1:end-1].+uInternal[ceil(Int64, ϕ*dims[1]):dims[1], 2:end], dims=2)
    end
    tmp = tmp[:,1]
    dν = 1.0/(dims[1]-1)
    M̃ϕ = 0.5*dν*sum(tmp[1:end-1].+tmp[2:end]) 

    return M̃ϕ
end

function M̃ϕ(u, W, dims, dν, hᵥ, ϕ; thresh="floor")

    uInternal = reshape(W*hᵥ*u, dims...)
    if thresh == "floor"
        M̃ϕ = sum(selectdim(uInternal, 1, floor(Int64, ϕ*dims[1]):dims[1]))
    else
        M̃ϕ = sum(selectdim(uInternal, 1, ceil(Int64, ϕ*dims[1]):dims[1]))
    end

    # uInternal = reshape(hᵥ*u, dims...)
    # dx = sqrt(π)/(dims[2]-1)
    # yMax = sqrt(π)
    # tmp = 0.5*yMax*dx.*sum(uInternal[ceil(Int64, ϕ*dims[1]):dims[1], 1:end-1].+uInternal[ceil(Int64, ϕ*dims[1]):dims[1], 2:end], dims=2)
    # tmp = tmp[:,1]
    # dν = 1.0/(dims[1]-1)
    # M̃ϕ = 0.5*dν*sum(tmp[1:end-1].+tmp[2:end]) 
    return M̃ϕ
end

# Dimensional bulk functional mass integrated over space and polymerisation 
# \ref{eq:Mstar}
function Mstarϕ(u, W, dims, dν, hᵥ, α_C, 𝒞, ϕ)
    return M̃ϕ(u, W, dims, dν, hᵥ, ϕ)*α_C*𝒞/(π*(1+α_C))
end

# function t̃_To_tStar(T̃ᵣ, N, ℰ, Ω, Ωperp, C_b, S_b, k₁, k₂, k₃, k_Ca, k_Cd, k_Sa, k_Sd)
function t̃_To_tStar(t̃, N, K₂, K₃, k₁, E₀, σ)
    return t̃*(N^2)*(K₂+σ*K₃)/(k₁*E₀)
end

# Where τ = t̃-t̃₀
function M̃ϕAnalytic(ϕ, ν₀, τ, α_C, β, K₂, K₄) 
    # Uniform thickness =>
    Ẽ = K₂/(1+K₂)
    a = erf( ((ϕ-ν₀)*(1+α_C) - Ẽ*β*τ)/sqrt( 4*Ẽ*K₂*K₄*(1+α_C)*τ) )
    return 0.5*π*(1.0-a)
end

function 𝒫star₅₀Analytic(h₀, h_C, h_S, k₁, k₂, k₃, k₄, Ωperp, 𝒮, 𝒞, ℰ, N, ϕ) 
    u = h₀/h_C
    λ = h_C/h_S
    ζ = (2*k₂*Ωperp)/(k₃*𝒮)
    γ = (2*k₂*Ωperp)/(k₁*𝒞)
    Δ = 2*k₂*k₄*Ωperp/(k₁*k₃*𝒮)
    F = (u*(1-Δ*(1+λ*u))) / ( (1+u)*(1 + ζ*(1+λ*u))*(1 + u + (1/γ))) 
    return ((k₁*𝒞*ℰ)/(4.0*Ωperp*N*ϕ))*F
end

# Eq 50
function 𝒫star₅₀Numeric(N, k₁, k₂, k₃, 𝒞, ℰ, 𝒮, h₀, k_Ca, k_Cd, k_Sa, k_Sd, Ωperp, T̃ᵣ₅₀)
    h_C = 2*k_Ca/k_Cd
    h_S = 2*k_Sa/k_Sd
    ζ = (2*k₂*Ωperp)/(k₃*𝒮)
    a = ((k₁*𝒞 )^2)*ℰ/(k₃*𝒮)
    b = 1/(4.0*Ωperp*N^2)
    c = ((h₀/h_C)*(1.0+(h₀/h_S)))/(((1+(h₀/h_C))^2)*(1+ζ*(1.0+(h₀/h_S))))
    d = 1.0/T̃ᵣ₅₀
    return a*b*c*d
end

function homogeneousWidthC(ν̃, K₂, K₄, α_C, β, t̃)
    Ẽ = K₂/(1+K₂)
    M = 1.0
    ξ = ν̃ - Ẽ*β*t̃/(1+α_C)
    D = Ẽ*K₂*K₄/(1+α_C)
    return (M/sqrt(4.0*π*D*t̃))*exp(-ξ^2/(4.0*D*t̃))
end



function T̃ᵣ₅₀Analytic(𝒞, 𝒮, ϕ, h₀, h_C, h_S, k₁, k₂, k₃, k₄, Ωperp, N, ν₀, t̃₀) 

    k₁*𝒞*ϕ

    u = h₀/h_C
    λ = h_C/h_S
    # ζ = (2*k₂*Ωperp)/(k₃*𝒮)
    γ = (2*k₂*Ωperp)/(k₁*𝒞)
    Δ = 2*k₂*k₄*Ωperp/(k₁*k₃*𝒮)

    a = k₁*𝒞*(ϕ-ν₀)
    b = (1 + γ*(1+u))*(1+λ*u)

    c = k₃*𝒮*γ*N 
    d = (1+u)*(1-Δ*(1+λ*u))

    return (t̃₀ + a*b/(c*d))
end


function T̃ᵣ₅₀Analytic2(K₂, ϕ, α_C, β, ν₀, t̃₀) 
    Ẽ = K₂/(1+K₂)
    return (t̃₀+(ϕ-ν₀)*(1+α_C)/(Ẽ*β))
end



export t̃_To_tStar
export M̃
export M̃ϕ
export M̃ϕ2
export Mstarϕ
export P_star
export homogeneousWidthC
export 𝒫star₅₀Analytic
export 𝒫star₅₀Numeric
export M̃ϕAnalytic
export T̃ᵣ₅₀Analytic
export T̃ᵣ₅₀Analytic2

end

# function 𝓟starUniform(𝒞, ℰ, 𝒮, ϕ, N, k₁, k₂, K₃, K₄, Ωperp, h₀, h_C, h_S)
#     Δ = k₁*𝒞/(2.0*k₂*Ωperp)
#     return (𝒞/(4*ϕ))*(k₁*ℰ/(Ωperp*N))* (h₀/((h₀+h_C)*(h₀+h_C*(1+Δ)))) * ((𝒮*Δ*K₃/𝒞 - K₄ - K₄*h₀/h_S)/(1 + 𝒮*Δ*K₃/𝒞 + h₀/h_S) )
# end
