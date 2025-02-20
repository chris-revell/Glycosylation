
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
using JLD2
using LinearAlgebra

@from "$(srcdir("Glycosylation.jl"))" using Glycosylation
@from "$(srcdir("UsefulFunctions.jl"))" using UsefulFunctions
@from "$(srcdir("DerivedParameters.jl"))" using DerivedParameters
@from "$(srcdir("CisternaWidth.jl"))" using CisternaWidth

#%%

subFolder = "Figure3"
folderName = "25-02-17-11-45-18_K₂=0.3_K₄=1.0_T̃ᵣ=0.385_differencing=centre_nSpatialDims=1_α_C=5.0_β=70.0_𝒟=204.0"

# # thicknessProfile = "Gaussian"
# differencing = "centre"
# solver = SSPRK432()
# nOutputs = 100
# # σGRF = 0.2
# σGaussian = 0.20

# nSpatialDims = 1
# Ngrid = 401
# dims = fill(Ngrid, nSpatialDims+1)

# include(projectdir("notebooks", "paramsRaw.jl"))

#%%

data1 = load(datadir("sims", subFolder, folderName, "solutionHVariation.jld2"))
@unpack sol1, p1, rawParams1 = data1
mat_h1 = reshape([p1.hᵥ[i,i] for i=1:prod(p1.dims)], p1.dims...)

data2 = load(datadir("sims", subFolder, folderName, "solutionFVariation.jld2"))
@unpack sol2, p2, rawParams2 = data2
mat_h2 = reshape([p2.hᵥ[i,i] for i=1:prod(p2.dims)], p2.dims...)

#%%

derivedParams = derivedParameters(rawParams1.Ω, rawParams1.Ωperp, rawParams1.N, rawParams1.k_Cd, rawParams1.k_Ca, rawParams1.k_Sd, rawParams1.k_Sa, rawParams1.k₁, rawParams1.k₂, rawParams1.k₃, rawParams1.k₄, rawParams1.𝒞, rawParams1.𝒮, rawParams1.ℰ, rawParams1.D_C, rawParams1.D_S, rawParams1.Tᵣstar; checks=true)
@unpack L₀, E₀, C_b, S_b, δ_C, δ_S, α_C, α_S, C₀, S₀, Tᵣ, T̃ᵣ, K₂, K₃, K₄, σ, ϵ, 𝒟, β, h_C, h_S, λ, ζ, γ, Δ, F = derivedParams

#%%

outLength = min(length(sol1.t), length(sol2.t))-5
frames = collect(1:outLength÷2-1:outLength)
fig = Figure(size=(1000,1000), fontsize=18)

νs = collect(range(0.0, 1.0, p1.dims[1]))
xMax = sqrt(π)
xs = collect(range(0.0, xMax, p1.dims[2]))

letterArray = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)", "(i)", "(j)", "(k)", "(l)", "(m)", "(n)", "(o)", "(p)"]

axesVec = [Axis(fig[1,1])]
lines!(axesVec[1], mat_h1[1,:], xs)
xlims!(axesVec[1], (0.0, 1.2*maximum(mat_h1[1,:])))
ylims!(axesVec[1], (0.0, xMax))
# axesVec[end].yticks = (0.0:sqrt(π):sqrt(π), [L"0.0", L"\sqrt{\pi}"])
axesVec[end].yticks = (0.0:sqrt(π):sqrt(π), [L"0.0", L"\sqrt{\pi}"])
axesVec[end].xticks = (0.0:1.0:1.0, [L"0.0", L"1.0"])
text!(axesVec[end], Point{2,Float64}(0.95*1.2*maximum(mat_h1[1,:]),1.6), text=popfirst!(letterArray), color=:black, align=(:right, :bottom)) 
for x=2:4
    uInternal = reshape(sol1.u[frames[x-1]], p1.dims...)
    push!(axesVec, Axis(fig[1,x]))
    # tString = @sprintf("%.3f", sol2.t[frames[x-1]])
    # Label(fig[1,x, Top()], L"\tilde{t} = %$tString")
    # Label(fig[1,x,Bottom()], popfirst!(letterArray))
    heatmap!(axesVec[end], νs, xs, uInternal, colormap=:batlow)
    text!(axesVec[end], Point{2,Float64}(0.95,1.6), text=popfirst!(letterArray), color=:white, align=(:right, :bottom)) 
end

for x=2:4  
    push!(axesVec, Axis(fig[2,x]))
    # Label(fig[2,x,Bottom()], popfirst!(letterArray))
    M = M̃(sol1.u[frames[x-1]], p1.W, p1.dims, p1.dν, p1.hᵥ)[:,1]
    lines!(axesVec[end], νs, M)
    text!(axesVec[end], Point{2,Float64}(0.95,(1.6/sqrt(π))*40.0), text=popfirst!(letterArray), color=:black, align=(:right, :bottom)) 
end

push!(axesVec, Axis(fig[3,1]))
# Label(fig[6,1,Top()], popfirst!(letterArray))
lines!(axesVec[end], p2.matFₑ, xs)
xlims!(axesVec[end], (0.0, 1.2*maximum(p2.matFₑ)))
ylims!(axesVec[end], (0.0, xMax))
axesVec[end].yticks = (0.0:sqrt(π):sqrt(π), [L"0.0", L"\sqrt{\pi}"])
axesVec[end].xticks = (0.0:1.0:1.0, [L"0.0", L"1.0"])
text!(axesVec[end], Point{2,Float64}(0.95*1.2*maximum(mat_h1[1,:]),1.6), text=popfirst!(letterArray), color=:black, align=(:right, :bottom)) 
for x=2:4
    uInternal = reshape(sol2.u[frames[x-1]], p2.dims...)
    push!(axesVec, Axis(fig[3,x]))
    # Label(fig[6,x,Top()], popfirst!(letterArray))
    heatmap!(axesVec[end], νs, xs, uInternal, colormap=:batlow)
    text!(axesVec[end], Point{2,Float64}(0.95,1.6), text=popfirst!(letterArray), color=:white, align=(:right, :bottom)) 
end

for x=2:4
    push!(axesVec, Axis(fig[4,x]))
    # Label(fig[8,x,Top()], popfirst!(letterArray))
    M = M̃(sol2.u[frames[x-1]], p2.W, p2.dims, p2.dν, p2.hᵥ)[:,1]
    lines!(axesVec[end], νs, M)
    text!(axesVec[end], Point{2,Float64}(0.95,(1.6/sqrt(π))*40.0), text=popfirst!(letterArray), color=:black, align=(:right, :bottom)) 
end

# axesVec[1].xlabel = L"h"
# axesVec[1].ylabel = L"x"
Label(fig[1,1,Left()], L"x")
Label(fig[1,1,BottomRight()], L"h")
axesVec[1].xgridvisible = false
for ax in axesVec[2:4]
    # ax.xlabel = L"\nu"
    # ax.ylabel = L"x"
    ax.xticks = (0.0:1.0:1.0, [L"0.0", L"1.0"])
    ax.yticks = (0.0:sqrt(π):sqrt(π), [L"0.0", L"\sqrt{\pi}"])    
    ax.xticklabelsvisible = false
    ax.yticklabelsvisible = false
    xlims!(ax, (0.0,1.0))
    ylims!(ax, (0.0,xMax))
end
# axesVec[5].ylabel = L"\tilde{M}"
axesVec[end].xgridvisible = false
Label(fig[2,2,Left()], L"\tilde{M}")
Label(fig[2,2,Bottom()], L"\nu")
Label(fig[2,3,Bottom()], L"\nu")
Label(fig[2,4,Bottom()], L"\nu")

for ax in axesVec[5:7]
    # ax.xlabel = L"\nu"
    # ax.ylabel = L"\tilde{M}"
    xlims!(ax, (0.0,1.0))
    ylims!(ax, (0.0,40.0))
    ax.xticks = (0.0:1.0:1.0, [L"0.0", L"1.0"])
    ax.yticks = (0.0:40.0:40.0, [L"0.0", L"40.0"])
    # ax.xticklabelsvisible = false
    ax.yticklabelsvisible = false
    ax.xgridvisible = false
    ax.ygridvisible = false
end
axesVec[5].yticklabelsvisible = true

# axesVec[8].xlabel = L"F_e"
# axesVec[8].ylabel = L"x"
Label(fig[3,1,Left()], L"x")
Label(fig[3,1,BottomRight()], L"F_e")
axesVec[8].xgridvisible = false
for ax in axesVec[9:11]
    # ax.xlabel = L"\nu"
    # ax.ylabel = L"x"
    ax.xticks = (0.0:1.0:1.0, [L"0.0", L"1.0"])
    ax.yticks = (0.0:sqrt(π):sqrt(π), [L"0.0", L"\sqrt{\pi}"])    
    ax.xticklabelsvisible = false
    ax.yticklabelsvisible = false
    xlims!(ax, (0.0,1.0))
    ylims!(ax, (0.0,xMax))
end
# axesVec[12].ylabel = L"\tilde{M}"
Label(fig[4,2,Left()], L"\tilde{M}")
Label(fig[4,2,Bottom()], L"\nu")
Label(fig[4,3,Bottom()], L"\nu")
Label(fig[4,4,Bottom()], L"\nu")
for ax in axesVec[12:end]
    # ax.xlabel = L"\nu"
    # ax.ylabel = L"\tilde{M}"
    xlims!(ax, (0.0,1.0))
    ylims!(ax, (0.0,40.0))
    ax.xticks = (0.0:1.0:1.0, [L"0.0", L"1.0"])
    ax.yticks = (0.0:40.0:40.0, [L"0.0", L"40.0"])
    # ax.xticklabelsvisible = false
    ax.yticklabelsvisible = false
end
axesVec[12].yticklabelsvisible = true

colsize!(fig.layout, 1, Aspect(1, 1.0))
colsize!(fig.layout, 2, Aspect(1, 1.0))
colsize!(fig.layout, 3, Aspect(1, 1.0))
colsize!(fig.layout, 4, Aspect(1, 1.0))
resize_to_layout!(fig)

display(fig)
save(datadir("sims",subFolder,folderName,"Figure3.png"), fig)

@show sol1.t[frames]
@show sol2.t[frames]