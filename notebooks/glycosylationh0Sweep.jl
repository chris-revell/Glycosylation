#%%
# flux_νₑ = (diffusive_flux_ν + advective_flux_ν)
# flux_νₑ = K₂*K₄.*Pν*∇ₑ*cᵥ - β*Pν*cₑ    where cᵥ is concentration over vertices, cₑ is concentration over edges 
# cₑ = Aᵤₚ*cᵥ
# flux_νₑ = (K₂*K₄.*Pν*∇ₑ - β*Pν*Aᵤₚ)*cᵥ
# flux_xyₑ = Dₑ*hₑ*diffusive_flux_xy
# flux_xyₑ = Dₑ*hₑ*Pxy*∇ₑ*cᵥ
# ċ = aE∇⋅flux_νₑ + a∇⋅flux_xyₑ
# ċ = a*E*∇⋅(K₂*K₄.*Pν*∇ₑ*cᵥ - β*Pν*Aᵤₚ*cᵥ) + a∇⋅(Dₑ*hₑ*Pxy*∇ₑ*cᵥ)
# Dₑ constant over edges 
# ċ = a*(E*∇⋅(K₂*K₄.*Pν*∇ₑ - β*Pν*Aᵤₚ) + 𝓓.*∇⋅(hₑ*Pxy*∇ₑ))*cᵥ

# L = -W⁻¹*Aᵀ*𝓓*l⁻¹*A .+ W⁻¹*Aᵀ*V*Aᵤₚ # Express model as a matrix operator 


# Cνν = W⁻¹*Aᵀ*Pν*l⁻¹*A
# Cν = Aᵀ*l⁻¹*Pν*Aᵤₚ
# flux_νₑ = (diffusive_flux_ν + advective_flux_ν)
# flux_νₑ = K₂*K₄.*Pν*∇ₑ*cᵥ - β*Pν*cₑ    where cᵥ is concentration over vertices, cₑ is concentration over edges 
# cₑ = Aᵤₚ*cᵥ
# flux_νₑ = (K₂*K₄.*Pν*∇ₑ - β*Pν*Aᵤₚ)*cᵥ
# flux_xyₑ = Dₑ*hₑ*diffusive_flux_xy
# flux_xyₑ = Dₑ*hₑ*Pxy*∇ₑ*cᵥ
# ċ = aE∇⋅flux_νₑ + a∇⋅flux_xyₑ
# ċ = a*E*∇⋅(K₂*K₄.*Pν*∇ₑ*cᵥ - β*Pν*Aᵤₚ*cᵥ) + a∇⋅(Dₑ*hₑ*Pxy*∇ₑ*cᵥ)
# Dₑ constant over edges 
# ċ = a*(E*∇⋅(K₂*K₄.*Pν*∇ₑ - β*Pν*Aᵤₚ) + 𝓓.*∇⋅(hₑ*Pxy*∇ₑ))*cᵥ


#

using OrdinaryDiffEq
using SparseArrays
using UnPack
using CairoMakie 
using FromFile
using DrWatson
using Printf
using SciMLOperators
using Dates
using InvertedIndices

@from "$(srcdir("Glycosylation.jl"))" using Glycosylation
@from "$(srcdir("Visualise.jl"))" using Visualise
@from "$(srcdir("UsefulFunctions.jl"))" using UsefulFunctions
@from "$(srcdir("MakeWeightMatrices.jl"))" using MakeWeightMatrices
@from "$(srcdir("DerivedParameters.jl"))" using DerivedParameters

#%%

include(projectdir("notebooks","paramsRaw.jl"))
h_C = 2*k_Ca/k_Cd
h_S = 2*k_Sa/k_Sd
hMax = h_C*100
hMin = h_C/10
h₀s = collect(hMin:5*hMin:hMax)
Ωs    = h₀s.*Ωperp      # Dimensional lumen volume 

#%%

# include(projectdir("notebooks","paramsDerived.jl"))

#%%

# h₀s = collect(0.01:0.01:1.0)
# Ωs    = h₀s.*Ωperp      # Dimensional lumen volume 

# Tᵣstar = (2\vert\Omega_\perp\vert N^2 \tilde{T}_r)/𝓔
#     \frac{2k_{Sa}\vert\Omega_\perp\vert+k_{Sd}\vert\Omega\vert}{2k_{Ca}\vert\Omega_\perp\vert+k_{Cd}\vert\Omega\vert}
#     \frac{k_1 k_{Ca} C_b \vert\Omega\vert}{k_1 k_2(2k_{Sa}\vert\Omega_\perp\vert+k_{Sd}\vert\Omega\vert)+k_3 k_{Sa} S_b\vert\Omega\vert}.




# Δ = k₁*𝓒/(2.0*k₂*Ωperp)

dims[2] = dims[2]÷100

xMax = π^(1/nSpatialDims)
xs   = collect(range(0.0, xMax, dims[2]))
dx   = xs[2]-xs[1]
if nSpatialDims > 1 
    yMax = xMax
    ys   = collect(range(0.0, yMax, dims[3]))
    dy   = ys[2]-ys[1]
end
νMax = 1.0
νs   = collect(range(0.0, νMax, dims[1])) # Positions of discretised vertices in polymerisation space 
dν   = νs[2]-νs[1]
nSpatialDims == 1 ? spacing  = [dν, dx] : spacing  = [dν, dx, dy]
W = vertexVolumeWeightsMatrix(dims, spacing)

PstarsAnalytic = []
PstarsSim = []
MstarsSim = []
sols = []
for i=1:length(h₀s)
    @show h₀s[i]
    derivedParams = derivedParameters(Ωs[i], Ωperp, N, k_Cd, k_Ca, k_Sd, k_Sa, k₁, k₂, k₃, k₄, 𝓒, 𝓢, 𝓔, D_C, D_S, Tᵣstar; checks=false)
    @unpack L₀, E₀, h₀, C_b, S_b, δ_C, δ_S, α_C, α_S, C₀, S₀, Tᵣ, K₂, K₃, K₄, σ, ϵ, 𝓓, β, λ = derivedParams

    # @show (𝓔*β*Tᵣ < 1+α_C)

    sol, mat_h = glycosylationAnyD(dims, K₂, K₄, Tᵣ, α_C, 𝓓, β, thickness="uniform", differencing="upstream", solver=SSPRK432())

    hᵥ = spdiagm(ones(prod(dims)))
    M_stars = Float64[]
    for u in sol.u
        uInternal = reshape(W*hᵥ*u, dims...)
        M̃ = sum(uInternal, dims=(2:length(dims)))
        Mϕ = sum(M̃[ceil(Int, ϕ*dims[1]) : dims[1]])
        push!(M_stars, Mϕ/sum(M̃))
    end
    push!(sols, sol)
    push!(PstarsSim, P_star(sol.u[end], W, dims, dν, hᵥ, α_C, C_b, Ωs[i], ϕ, Ωperp, k₁, 𝓔, Tᵣ))
    # push!(PstarsSim, M_stars[end]/sol.t[end])
    push!(MstarsSim, M_stars[end])
    push!(PstarsAnalytic, 𝓟starUniform(𝓒, 𝓔, 𝓢, ϕ, N, k₁, k₂, K₃, K₄, Ωperp, h₀s[i], h_C, h_S))
end


# PstarsAnalytic = []
# PstarsSim = []
# j = 5
# for (i,sol) in enumerate(sols)
#     derivedParams = derivedParameters(Ωs[i], Ωperp, N, k_Cd, k_Ca, k_Sd, k_Sa, k₁, k₂, k₃, k₄, 𝓒, 𝓢, 𝓔, D_C, D_S, Tᵣstar; checks=false)
#     @unpack L₀, E₀, h₀, C_b, S_b, δ_C, δ_S, α_C, α_S, C₀, S₀, Tᵣ, K₂, K₃, K₄, σ, ϵ, 𝓓, β, λ = derivedParams
#     # push!(PstarsSim, P_star(sol.u[j], W, dims, dν, spdiagm(ones(prod(dims))), α_C, C_b, Ωs[i], ϕ, sol.t[j]))
#     # push!(MstarsSim, M_star_ϕ(sol.u[end], W, dims, dν, spdiagm(ones(prod(dims))), α_C, C_b, Ωs[i], ϕ))
#     push!(PstarsAnalytic, 𝓟starUniform(𝓒, 𝓔, 𝓢, ϕ, N, k₁, K₃, K₄, Ωperp, h₀s[i], h_C, h_S, Δ))
# end

#%%

linesVec = []
labelsVec = []
fig = Figure(size=(500,500))
ax1 = Axis(fig[1,1])
push!(linesVec, lines!(ax1, h₀s, MstarsSim, color=:blue))
push!(labelsVec, "Numerical")
ylims!(ax1, (0.0,maximum(MstarsSim)))
xlims!(ax1, (0.0,maximum(h₀s)))
ax1.xlabel = "h₀"
ax1.ylabel = L"M^*"

ax2 = Axis(fig[2,1])
push!(linesVec, lines!(ax2, h₀s, PstarsSim, color=:blue))
push!(labelsVec, "Numerical")
ylims!(ax2, (0.0,maximum(PstarsSim)))
xlims!(ax2, (0.0,maximum(h₀s)))
ax2.xlabel = "h₀"
ax2.ylabel = L"𝓟^*"

ax3 = Axis(fig[3,1])
ax3.xlabel = "h₀"
ax3.ylabel = L"𝓟^*"
ylims!(ax3, (0.0,maximum(PstarsAnalytic)))
xlims!(ax3, (0.0,maximum(h₀s)))
push!(linesVec, lines!(ax3, h₀s, PstarsAnalytic, color=:red))
push!(labelsVec, "Analytic")

# push!(linesVec, vlines!(ax1, h_C, color=:green))
# push!(labelsVec, L"h_C")
# push!(linesVec, vlines!(ax2, h_C, color=:green))
# push!(labelsVec, L"h_C")
# push!(linesVec, vlines!(ax1, h_S, color=:orange))
# push!(labelsVec, L"h_S")
# push!(linesVec, vlines!(ax2, h_S, color=:orange))
# push!(labelsVec, L"h_S")

# Legend(fig[:,2], linesVec, labelsVec)

display(fig)

paramsName = @savename k_Cd k_Ca k_Sd k_Sa k₁ k₂ k₃ k₄ 𝓒 𝓢 𝓔 D_C D_S Tᵣstar ϕ
# folderName = "$(Dates.format(Dates.now(),"yy-mm-dd-HH-MM-SS"))_$(paramsName)"
save("$(Dates.format(Dates.now(),"yy-mm-dd-HH-MM-SS"))_$(paramsName)_simulationPvsh.png",fig)



