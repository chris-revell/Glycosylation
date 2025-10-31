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

differencing = "centre"
nSpatialDims = 1
T̃ᵣ = 30.0
K₂ = 1.0
K₄ = 0.0001
α_C = 1.0
𝓓 = 1.0
β = 0.1
Ngrid = 401
# dims  = fill(Ngrid, nSpatialDims+1)
dims = [Ngrid,2]

#%%

# sol = glycosylation(dims, K₂, K₄, T̃ᵣ, α_C, 𝓓, β, thickness="uniform", differencing=differencing) 
solCentre = glycosylation(dims, K₂, K₄, T̃ᵣ, α_C, 𝓓, β, thickness="uniform", differencing="centre", solver=SSPRK432())#NDBLSRK124()) 
solUpstream = glycosylation(dims, K₂, K₄, T̃ᵣ, α_C, 𝓓, β, thickness="uniform", differencing="upstream", solver=SSPRK432())#NDBLSRK124()) 
println("finished sim")

#%%

# Create directory for run data labelled with current time.
paramsName = @savename nSpatialDims K₂ K₄ α_C β 𝓓 T̃ᵣ differencing
folderName = "$(Dates.format(Dates.now(),"yy-mm-dd-HH-MM-SS"))_$(paramsName)"
subFolder = "analyticNumericFit"
mkpath(datadir("sims",subFolder,folderName))

#%%

νs   = collect(range(0.0, 1.0, dims[1])) 

# midpoint = length(sol.u)÷2
# C_peak, ind_peak = findmax(reshape(sol.u[midpoint], dims...)[:,1])
# ν_peak = νs[ind_peak]
# Ẽ = K₂/(1+K₂)
# D = Ẽ*K₂*K₄/(1+α_C)
# t₀ = sol.t[midpoint] - 1/(4.0*π*D*C_peak^2)
# ν₀ = ν_peak - Ẽ*β*(sol.t[midpoint]-t₀)/(1+α_C)

# νsOffset = νs.-ν₀
# tsOffset = sol.t.-t₀

#%%

# fig = Figure(size=(1000,750), fontsize=32)
# ax = CairoMakie.Axis(fig[1, 1])
# ylims!(ax, (0.0, 20.0))
# xlims!(ax, (0.0, 1.0))
# allLines = []
# for (c,i) in enumerate([1, 251, 500])
#     uInternal = reshape(sol.u[i], dims...)
#     push!(allLines, lines!(ax, νs, uInternal[:,1], linestyle=:solid, color=(:blue, 1.0), linewidth=4))
#     push!(allLines, lines!(ax, νs, homogeneousWidthC.(νsOffset, K₂, K₄, α_C, β, tsOffset[i]), linestyle=:dash, color=(:red,1.0), linewidth=4))
# end
# Legend(fig[1,2], allLines[1:2], ["Numeric", "Analytic"])
# ax.xlabel = L"\nu"
# ax.ylabel = L"C"
# save(datadir("sims", subFolder, folderName, "analyticComparison.png"), fig)

# #%%

# fig = Figure(size=(1000,1000))
# ax = CairoMakie.Axis(fig[1, 1], aspect=1)
# ax.xlabel = "ν"
# ax.ylabel = "C"
# analyticLine = Observable(zeros(dims[1]))
# numericLine = Observable(zeros(dims[1]))
# l1 = lines!(ax, νs, analyticLine, color=:red)
# l2 = lines!(ax, νs, numericLine, color=:blue)
# Legend(fig[1,2], [l1, l2], ["Analytic", "Numeric"])
# ylims!(ax, (-20.0, 20.0))
# xlims!(ax, (0.0, 1.0))
# analyticVals = homogeneousWidthC.(νsOffset, K₂, K₄, α_C, β, tsOffset[1])
# record(fig, datadir("sims",subFolder, folderName, "analyticCs.mp4"), 1:length(sol.t); framerate=50) do i
#     analyticVals .= homogeneousWidthC.(νsOffset, K₂, K₄, α_C, β, tsOffset[i])
#     analyticLine[] .= analyticVals
#     uInternal = reshape(sol.u[i], dims...)
#     numericLine[] .= uInternal[:,dims[2]÷2]
#     analyticLine[] = analyticLine[]
#     numericLine[] = numericLine[]
#     # if i in [1, 251, 500]
#     #     save(datadir("sims",subFolder, folderName, "analyticCs$i.png"), fig)
#     # end
# end

#%%

fig = Figure(size=(1000,1000))
ax = CairoMakie.Axis(fig[1, 1], aspect=1)
ax.xlabel = "ν"
ax.ylabel = "C"
analyticLine = Observable(zeros(dims[1]))
numericLine = Observable(zeros(dims[1]))
l1 = lines!(ax, νs, analyticLine, color=(:red, 0.5), linewidth=4, linestyle=:dot)
l2 = lines!(ax, νs, numericLine, color=(:blue, 0.5), linewidth=4)#, linestyle=:dot)
Legend(fig[1,2], [l1, l2], ["Upstream differencing", "Central differencing"])
ylims!(ax, (-20.0, 20.0))
xlims!(ax, (0.0, 1.0))
# analyticVals = homogeneousWidthC.(νsOffset, K₂, K₄, α_C, β, tsOffset[1])
record(fig, datadir("sims",subFolder, folderName, "analyticCs.mp4"), 1:length(solUpstream.t); framerate=50) do i
    uInternal = reshape(solUpstream.u[i], dims...)
    analyticLine[] .= uInternal[:,1]
    uInternal = reshape(solCentre.u[i], dims...)
    numericLine[] .= uInternal[:,1]
    analyticLine[] = analyticLine[]
    numericLine[] = numericLine[]
    # if i in [1, 251, 500]
    #     save(datadir("sims",subFolder, folderName, "analyticCs$i.png"), fig)
    # end
end
