
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
using XLSX
using DataFrames
using Interpolations
using Statistics
using JLD2

@from "$(srcdir("Glycosylation.jl"))" using Glycosylation
@from "$(srcdir("Visualise.jl"))" using Visualise
@from "$(srcdir("UsefulFunctions.jl"))" using UsefulFunctions
@from "$(srcdir("MakeWeightMatrices.jl"))" using MakeWeightMatrices
@from "$(srcdir("DerivedParameters.jl"))" using DerivedParameters
@from "$(srcdir("CisternaWidth.jl"))" using CisternaWidth

#%%

thicknessProfile = "Gaussian"
differencing = "centre"
solver = SSPRK432()
nOutputs = 100

σGaussian = 0.20

nSpatialDims = 1
Ngrid = 401
dims = fill(Ngrid, nSpatialDims+1)

#%%

h₀ = 0.1
Ωperp = 10000    # Dimensional lumen footprint area
Ω     = h₀*Ωperp      # Dimensional lumen volume 
N     = 100     # Maximum polymer length 
k_Cd  = 1.0 # Complex desorption rate
k_Ca  = 0.01 # Complex adsorption rate
k_Sd  = 1.0 # Substrate desorption rate
k_Sa  = 0.01 # Substrate adsorption rate
k₁    = 1.0   # Complex formation forward reaction rate 
k₂    = 0.1   # Complex dissociation reverse reaction rate 
k₃    = 0.1   # Product formation
k₄    = 0.1  # Product dissociation 
𝓒     = 100000.0
𝓢     = 100000.0
𝓔     = 0.0001
D_C   = 0.0000001  # Monomer/polymer diffusivity
D_S   = 0.0000001  # Substrate diffusivity
Tᵣstar= 1000000000.0  # Release time
ϕ     = 0.5

xMax = π^(1/nSpatialDims)
xs   = collect(range(0.0, xMax, dims[2]))
# dx   = xs[2]-xs[1]
if nSpatialDims > 1 
    yMax = xMax
    ys   = collect(range(0.0, yMax, dims[3]))
    # dy   = ys[2]-ys[1]
end
νMax = 1.0
νs   = collect(range(0.0, νMax, dims[1]))
# dν   = νs[2]-νs[1]
# nSpatialDims == 1 ? spacing  = [dν, dx] : spacing  = [dν, dx, dy]
# W = vertexVolumeWeightsMatrix(dims, spacing)

#%%

derivedParams = derivedParameters(Ω, Ωperp, N, k_Cd, k_Ca, k_Sd, k_Sa, k₁, k₂, k₃, k₄, 𝓒, 𝓢, 𝓔, D_C, D_S, Tᵣstar; checks=true)
@unpack L₀, E₀, h₀, C_b, S_b, δ_C, δ_S, α_C, α_S, C₀, S₀, Tᵣ, K₂, K₃, K₄, σ, ϵ, 𝓓, β, λ = derivedParams

#%%


# Create directory for run data labelled with current time.
# paramsName = @savename nSpatialDims K₂ K₃ K₄ α_C δ_C σ N β 𝓓 Tᵣ h₀ Ωperp 𝓒
paramsName = @savename nSpatialDims K₂ K₄ α_C β 𝓓 Tᵣ differencing
folderName = "$(Dates.format(Dates.now(),"yy-mm-dd-HH-MM-SS"))_$(paramsName)"
# Create frames subdirectory to store system state at each output time
subFolder = "VaryingHandF"
mkpath(datadir("sims",subFolder,folderName))

#%%

sol1, p1 = glycosylationAnyD(dims, K₂, K₄, Tᵣ, α_C, 𝓓, β, thickness=thicknessProfile, differencing=differencing, solver=solver, nOutputs=nOutputs, σGaussian=σGaussian)
println("finished sim")

# @unpack hᵥ, matFₑ, W = p
mat_h = reshape([p1.hᵥ[i,i] for i=1:prod(dims)], dims...)

#%%

jldsave(datadir("sims",subFolder,folderName,"solutionHVariation.jld2"); sol1, p1)

#%%

# globalmin = minimum([minimum(u) for u in sol1.u[4:end]])
# globalmax = maximum([maximum(u) for u in sol1.u[4:end]])
# clims = (globalmin,globalmax)
# # frames = collect(4:48:100)
# frames = collect(1:40:100)
# fig = Figure()#size=(1000,500))
# axesVec = [Axis(fig[1,1])]
# lines!(axesVec[1], mat_h[1,:], xs)
# xlims!(axesVec[1], (0.0, 1.2*maximum(mat_h[1,:])))
# ylims!(axesVec[1], (0.0, xMax))
# for x=2:4
#     uInternal = reshape(sol1.u[frames[x-1]], dims...)
#     push!(axesVec, Axis(fig[1,x]))
#     axesVec[end].title = "t = $(@sprintf("%.2f", sol1.t[frames[x-1]]))"
#     # heatmap!(axesVec[end], νs, xs, uInternal, colorrange=clims, colormap=:batlow)
#     heatmap!(axesVec[end], νs, xs, uInternal, colormap=:batlow)
#     xlims!(axesVec[end], (0.0,1.0))
#     ylims!(axesVec[end], (0.0,xMax))
#     push!(axesVec, Axis(fig[2,x]))
#     M = M_tilde(sol1.u[frames[x-1]], p1.W, dims, p1.dν, p1.hᵥ)[:,1]
#     lines!(axesVec[end], νs, M)
#     xlims!(axesVec[end], (0.0,1.0))
# end
# axesVec[1].xlabel = L"h"
# axesVec[1].ylabel = "x"
# for ax in axesVec[2:2:end]
#     ax.xlabel = L"\nu"
#     ax.ylabel = L"x"
# end
# for ax in axesVec[3:2:end]
#     ax.xlabel = L"\nu"
#     ax.ylabel = L"M"
# end
# # display(fig)
# save(datadir("sims",subFolder,folderName,"2DhVariation.png"), fig)

#%%

sol2, p2 = glycosylationAnyD(dims, K₂, K₄, Tᵣ, α_C, 𝓓, β, thickness="uniform", fDist="Gaussian", differencing=differencing, solver=solver, nOutputs=nOutputs, σGaussian=σGaussian)
println("finished sim 2")

# @unpack hᵥ, matFₑ, W = p2
# mat_h2 = reshape([p2.hᵥ[i,i] for i=1:prod(dims)], dims...)

#%%

jldsave(datadir("sims",subFolder,folderName,"solutionFVariation.jld2"); sol2, p2)

#%%

# globalmin = minimum([minimum(u) for u in sol2.u[4:end]])
# globalmax = maximum([maximum(u) for u in sol2.u[4:end]])
# clims = (globalmin,globalmax)
# # frames = collect(4:48:100)
# frames = collect(1:40:100)
# fig = Figure()#size=(1000,500))
# axesVec = [Axis(fig[1,1])]
# lines!(axesVec[1], p2.matFₑ, xs)
# xlims!(axesVec[1], (0.0, 1.2*maximum(p2.matFₑ)))
# ylims!(axesVec[1], (0.0, xMax))
# axesVec[1].xlabel = L"F_e"
# axesVec[1].ylabel = L"x"
# for x=2:4
#     uInternal = reshape(sol2.u[frames[x-1]], dims...)
#     push!(axesVec, Axis(fig[1,x]))
#     axesVec[end].title = "t = $(@sprintf("%.2f", sol2.t[frames[x-1]]))"
#     # heatmap!(axesVec[end], νs, xs, uInternal, colorrange=clims, colormap=:batlow)
#     heatmap!(axesVec[end], νs, xs, uInternal, colormap=:batlow)
#     xlims!(axesVec[end], (0.0,1.0))
#     ylims!(axesVec[end], (0.0,xMax))
#     push!(axesVec, Axis(fig[2,x]))
#     M = M_tilde(sol2.u[frames[x-1]], p.W, dims, p2.dν, p2.hᵥ)[:,1]
#     lines!(axesVec[end], νs, M)
#     xlims!(axesVec[end], (0.0,1.0))
# end
# for ax in axesVec[2:2:end]
#     ax.xlabel = L"\nu"
#     ax.ylabel = L"x"
# end
# for ax in axesVec[3:2:end]
#     ax.xlabel = L"\nu"
#     ax.ylabel = L"\tilde{M}"
# end
# # display(fig)
# save(datadir("sims",subFolder,folderName,"2DFVariation.png"), fig)

#%%

# frames = collect(4:48:100)
frames = collect(1:40:100)
fig = Figure(size=(1000,1000))

# globalmin1 = minimum([minimum(u) for u in sol1.u[4:end]])
# globalmax1 = maximum([maximum(u) for u in sol1.u[4:end]])
# # # clims1 = (globalmin1,globalmax1)
axesVec = [Axis(fig[1,1])]
lines!(axesVec[1], mat_h[1,:], xs)
xlims!(axesVec[1], (0.0, 1.2*maximum(mat_h[1,:])))
ylims!(axesVec[1], (0.0, xMax))
for x=2:4
    uInternal = reshape(sol1.u[frames[x-1]], dims...)
    push!(axesVec, Axis(fig[1,x]))
    # axesVec[end].title = "t = $(@sprintf("%.2f", sol1.t[frames[x-1]]))"
    Label(fig[1,x, Top()], "t = $(@sprintf("%.2f", sol1.t[frames[x-1]]))")
    # heatmap!(axesVec[end], νs, xs, uInternal, colorrange=clims, colormap=:batlow)
    heatmap!(axesVec[end], νs, xs, uInternal, colormap=:batlow)
    xlims!(axesVec[end], (0.0,1.0))
    ylims!(axesVec[end], (0.0,xMax))
    push!(axesVec, Axis(fig[2,x]))
    M = M_tilde(sol1.u[frames[x-1]], p1.W, dims, p1.dν, p1.hᵥ)[:,1]
    lines!(axesVec[end], νs, M)
    xlims!(axesVec[end], (0.0,1.0))
end
axesVec[1].xlabel = L"h"
axesVec[1].ylabel = L"x"
for ax in axesVec[2:2:end]
    ax.xlabel = L"\nu"
    ax.ylabel = L"x"
end
for ax in axesVec[3:2:end]
    ax.xlabel = L"\nu"
    ax.ylabel = L"\tilde{M}"
end

# globalmin2 = minimum([minimum(u) for u in sol2.u[4:end]])
# globalmax2 = maximum([maximum(u) for u in sol2.u[4:end]])
# # # clims2 = (globalmin2,globalmax2)
push!(axesVec, Axis(fig[3,1]))
lines!(axesVec[end], p2.matFₑ, xs)
xlims!(axesVec[end], (0.0, 1.2*maximum(p2.matFₑ)))
ylims!(axesVec[end], (0.0, xMax))
axesVec[end].xlabel = L"F_e"
axesVec[end].ylabel = L"x"
for x=2:4
    uInternal = reshape(sol2.u[frames[x-1]], dims...)
    push!(axesVec, Axis(fig[3,x]))
    # axesVec[end].title = "t = $(@sprintf("%.2f", sol2.t[frames[x-1]]))"
    Label(fig[3,x, Top()], "t = $(@sprintf("%.2f", sol2.t[frames[x-1]]))")
    # heatmap!(axesVec[end], νs, xs, uInternal, colorrange=clims, colormap=:batlow)
    heatmap!(axesVec[end], νs, xs, uInternal, colormap=:batlow)
    xlims!(axesVec[end], (0.0,1.0))
    ylims!(axesVec[end], (0.0,xMax))
    push!(axesVec, Axis(fig[4,x]))
    M = M_tilde(sol2.u[frames[x-1]], p.W, dims, p2.dν, p2.hᵥ)[:,1]
    lines!(axesVec[end], νs, M)
    xlims!(axesVec[end], (0.0,1.0))
end
for ax in axesVec[9:2:end]
    ax.xlabel = L"\nu"
    ax.ylabel = L"x"
end
for ax in axesVec[10:2:end]
    ax.xlabel = L"\nu"
    ax.ylabel = L"\tilde{M}"
end

linkxaxes!(axesVec[2], axesVec[3:7])
linkxaxes!(axesVec[2], axesVec[9:end])
linkyaxes!(axesVec[1], axesVec[2:2:7])
linkyaxes!(axesVec[1], axesVec[8])
linkyaxes!(axesVec[1], axesVec[9:2:end])
linkyaxes!(axesVec[3], axesVec[5:2:7])
linkyaxes!(axesVec[10], axesVec[12:2:end])

display(fig)
save(datadir("sims",subFolder,folderName,"2DFandHVariation.png"), fig)
