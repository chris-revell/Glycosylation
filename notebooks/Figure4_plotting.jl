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
# @from "$(srcdir("CisternaWidth.jl"))" using CisternaWidth

subFolder = "Figure4"
folderName = "25-11-07-17-16-06_K₂=0.3_K₄=1.0_T̃ᵣ=3.85_differencing=centre_nSpatialDims=1_α_C=5.0_β=70.0_𝒟=204.0"

data0 = load(datadir("sims", subFolder, folderName, "solutionNoVariation.jld2"))
@unpack sol0, p0 = data0
mat_h0 = reshape([p0.hᵥ[i,i] for i=1:prod(p0.dims)], p0.dims...)

data1 = load(datadir("sims", subFolder, folderName, "solutionHVariation.jld2"))
@unpack sol1, p1 = data1
mat_h1 = reshape([p1.hᵥ[i,i] for i=1:prod(p1.dims)], p1.dims...)

data2 = load(datadir("sims", subFolder, folderName, "solutionFVariation.jld2"))
@unpack sol2, p2 = data2
mat_h2 = reshape([p2.hᵥ[i,i] for i=1:prod(p2.dims)], p2.dims...)

# derivedParams = derivedParameters(rawParams1.Ω, rawParams1.𝒜, rawParams1.N, rawParams1.k_Cd, rawParams1.k_Ca, rawParams1.k_Sd, rawParams1.k_Sa, rawParams1.k₁, rawParams1.k₂, rawParams1.k₃, rawParams1.k₄, rawParams1.𝒞, rawParams1.𝒮, rawParams1.ℰ, rawParams1.D_C, rawParams1.D_S, rawParams1.Tᵣstar; checks=true)
# @unpack L₀, E₀, C_b, S_b, δ_C, δ_S, α_C, α_S, C₀, S₀, Tᵣ, T̃ᵣ, K₂, K₃, K₄, σ, ϵ, 𝒟, β, h_C, h_S, λ, ζ, γ, Δ, F = derivedParams

outLength = minimum([length(sol0.t), length(sol1.t), length(sol2.t)])-15
frames = collect(1:outLength÷2-1:outLength)
fig = Figure(size=(1000,1000), fontsize=18)

νs = collect(range(0.0, 1.0, p1.dims[1]))
xMax = sqrt(π)
xs = collect(range(0.0, xMax, p1.dims[2]))

letterArray = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P"]

axesVec = [Axis(fig[1,1])]
lines!(axesVec[1], mat_h1[1,:], xs, linewidth=2)
xlims!(axesVec[1], (0.0, 1.2*maximum(mat_h1[1,:])))
ylims!(axesVec[1], (0.0, xMax))
# axesVec[end].yticks = (0.0:sqrt(π):sqrt(π), [L"0.0", L"\sqrt{\pi}"])
axesVec[end].yticks = (0.0:sqrt(π):sqrt(π), [L"0.0", L"\sqrt{\pi}"])
axesVec[end].xticks = (0.0:1.0:1.0, [L"0.0", L"1.0"])
text!(axesVec[end], Point{2,Float64}(0.95*1.2*maximum(mat_h1[1,:]),1.5), text=popfirst!(letterArray), color=:black, align=(:right, :bottom), fontsize=24) 
for x=2:4
    uInternal = reshape(sol1.u[frames[x-1]], p1.dims...)
    push!(axesVec, Axis(fig[1,x]))
    # tString = @sprintf("%.3f", sol2.t[frames[x-1]])
    # Label(fig[1,x, Top()], L"\tilde{t} = %$tString")
    # Label(fig[1,x,Bottom()], popfirst!(letterArray))
    heatmap!(axesVec[end], νs, xs, uInternal, colormap=:batlow)
    text!(axesVec[end], Point{2,Float64}(0.95,1.5), text=popfirst!(letterArray), color=:white, align=(:right, :bottom), fontsize=24) 
end

for x=2:4  
    push!(axesVec, Axis(fig[2,x]))
    # Label(fig[2,x,Bottom()], popfirst!(letterArray))
    M = M̃(sol1.u[frames[x-1]], p1.W, p1.dims, p1.dν, p1.hᵥ)[:,1]
    lines!(axesVec[end], νs, M, linewidth=2)
    @show sol1.t[frames[x-1]]
    Muniform = M̃(sol0.u[frames[x-1]], p0.W, p0.dims, p0.dν, p0.hᵥ)[:,1]
    lines!(axesVec[end], νs, Muniform, linewidth=2)
    @show sol0.t[frames[x-1]]
    text!(axesVec[end], Point{2,Float64}(0.95,(1.5/sqrt(π))*40.0), text=popfirst!(letterArray), color=:black, align=(:right, :bottom), fontsize=24)
end
vlines!(axesVec[end], 0.5, color=(:black,0.5))

push!(axesVec, Axis(fig[4,1]))
# Label(fig[6,1,Top()], popfirst!(letterArray))
lines!(axesVec[end], p2.matFₑ, xs, linewidth=2)
xlims!(axesVec[end], (0.0, 1.2*maximum(p2.matFₑ)))
ylims!(axesVec[end], (0.0, xMax))
axesVec[end].yticks = (0.0:sqrt(π):sqrt(π), [L"0.0", L"\sqrt{\pi}"])
axesVec[end].xticks = (0.0:1.0:1.0, [L"0.0", L"1.0"])
text!(axesVec[end], Point{2,Float64}(0.95*1.2*maximum(mat_h1[1,:]),1.5), text=popfirst!(letterArray), color=:black, align=(:right, :bottom), fontsize=24) 
for x=2:4
    uInternal = reshape(sol2.u[frames[x-1]], p2.dims...)
    push!(axesVec, Axis(fig[4,x]))
    # Label(fig[6,x,Top()], popfirst!(letterArray))
    heatmap!(axesVec[end], νs, xs, uInternal, colormap=:batlow)
    text!(axesVec[end], Point{2,Float64}(0.95,1.5), text=popfirst!(letterArray), color=:white, align=(:right, :bottom), fontsize=24) 
end

for x=2:4
    push!(axesVec, Axis(fig[5,x]))
    M = M̃(sol2.u[frames[x-1]], p2.W, p2.dims, p2.dν, p2.hᵥ)[:,1]
    lines!(axesVec[end], νs, M, linewidth=2)
    @show sol2.t[frames[x-1]]
    Muniform = M̃(sol0.u[frames[x-1]], p0.W, p0.dims, p0.dν, p0.hᵥ)[:,1]
    lines!(axesVec[end], νs, Muniform, linewidth=2)
    @show sol0.t[frames[x-1]]
    text!(axesVec[end], Point{2,Float64}(0.95,(1.5/sqrt(π))*40.0), text=popfirst!(letterArray), color=:black, align=(:right, :bottom), fontsize=24) 
end
vlines!(axesVec[end], 0.5, color=(:black,0.5))

Label(fig[1,1, Left()], L"x")
Label(fig[2,1, Top()], L"h")
axesVec[1].xgridvisible = false
for ax in axesVec[2:4]
    ax.xticks = (0.0:1.0:1.0, [L"0.0", L"1.0"])
    ax.yticks = (0.0:sqrt(π):sqrt(π), [L"0.0", L"\sqrt{\pi}"])    
    ax.xticklabelsvisible = false
    ax.yticklabelsvisible = false
    xlims!(ax, (0.0,1.0))
    ylims!(ax, (0.0,xMax))
end
axesVec[end].xgridvisible = false
Label(fig[2,2,Left()], L"\tilde{M}")
Label(fig[3,2,Top()], L"\nu")
Label(fig[3,3,Top()], L"\nu")
Label(fig[3,4,Top()], L"\nu")

for ax in axesVec[5:7]
    xlims!(ax, (0.0,1.0))
    ylims!(ax, (0.0,40.0))
    ax.xticks = (0.0:1.0:1.0, [L"0.0", L"1.0"])
    ax.yticks = (0.0:40.0:40.0, [L"0.0", L"40.0"])
    ax.yticklabelsvisible = false
    ax.xgridvisible = false
    ax.ygridvisible = false
end
axesVec[5].yticklabelsvisible = true
axesVec[7].xticks = (0.0:0.5:1.0, [L"0.0", L"\nu=\phi", L"1.0"])

Label(fig[4,1,Left()], L"x")
Label(fig[5,1,Top()], L"F_E")
axesVec[8].xgridvisible = false
for ax in axesVec[9:11]
    ax.xticks = (0.0:1.0:1.0, [L"0.0", L"1.0"])
    ax.yticks = (0.0:sqrt(π):sqrt(π), [L"0.0", L"\sqrt{\pi}"])    
    ax.xticklabelsvisible = false
    ax.yticklabelsvisible = false
    xlims!(ax, (0.0,1.0))
    ylims!(ax, (0.0,xMax))
end
Label(fig[5,2,Left()], L"\tilde{M}")
Label(fig[6,2,Top()], L"\nu")
Label(fig[6,3,Top()], L"\nu")
Label(fig[6,4,Top()], L"\nu")
for ax in axesVec[12:end]
    xlims!(ax, (0.0,1.0))
    ylims!(ax, (0.0,40.0))
    ax.xticks = (0.0:1.0:1.0, [L"0.0", L"1.0"])
    ax.yticks = (0.0:40.0:40.0, [L"0.0", L"40.0"])
    ax.yticklabelsvisible = false
end
axesVec[12].yticklabelsvisible = true
axesVec[end].xticks = (0.0:0.5:1.0, [L"0.0", L"\nu=\phi", L"1.0"])

colsize!(fig.layout, 1, Aspect(1, 1.0))
colsize!(fig.layout, 2, Aspect(1, 1.0))
colsize!(fig.layout, 3, Aspect(1, 1.0))
colsize!(fig.layout, 4, Aspect(1, 1.0))

rowsize!(fig.layout, 1, Relative(0.2475))
rowsize!(fig.layout, 2, Relative(0.2475))
rowsize!(fig.layout, 3, Relative(0.005))
rowsize!(fig.layout, 4, Relative(0.2475))
rowsize!(fig.layout, 5, Relative(0.2475))
rowsize!(fig.layout, 6, Relative(0.005))

rowgap!(fig.layout, 1, Relative(-0.003))
rowgap!(fig.layout, 2, Relative(-0.003))
rowgap!(fig.layout, 4, Relative(-0.003))
rowgap!(fig.layout, 5, Relative(-0.003))

resize_to_layout!(fig)

display(fig)
save(datadir("sims",subFolder,folderName,"Figure4.png"), fig)
save(datadir("sims",subFolder,folderName,"Figure4.pdf"), fig)

# @show sol1.t[frames]
# @show sol2.t[frames]

# include(projectdir("notebooks", "paramsRaw.jl"))
# derivedParams = derivedParameters(Ω, 𝒜, N, k_Cd, k_Ca, k_Sd, k_Sa, k₁, k₂, k₃, k₄, 𝒞, 𝒮, ℰ, D_C, D_S, Tᵣstar; checks=true)
# @unpack L₀, E₀, C_b, S_b, δ_C, δ_S, α_C, α_S, C₀, S₀, Tᵣ, T̃ᵣ, K₂, K₃, K₄, σ, ϵ, 𝒟, β, h_C, h_S, λ, ζ, γ, Δ, F = derivedParams
# ind50 = findfirst(x->M̃ϕ(x, p1.W, p1.dims, p1.dν, p1.hᵥ, 0.5) > 0.5*π, sol1.u)
# Tᵣ₅₀Star = sol1.t[ind50]*(N^2)*(K₂+σ*K₃)/(k₁*E₀)    
# 𝒫sim1 = Mstarϕ(sol1.u[ind50], p1.W, p1.dims, p1.dν, p1.hᵥ, α_C, 𝒞, 0.5)/Tᵣ₅₀Star
# @show 𝒫sim1

# ind50 = findfirst(x->M̃ϕ(x, p2.W, p2.dims, p2.dν, p2.hᵥ, 0.5) > 0.5*π, sol2.u)
# Tᵣ₅₀Star = sol2.t[ind50]*(N^2)*(K₂+σ*K₃)/(k₁*E₀)    
# 𝒫sim2 = Mstarϕ(sol2.u[ind50], p2.W, p2.dims, p2.dν, p2.hᵥ, α_C, 𝒞, 0.5)/Tᵣ₅₀Star
# @show 𝒫sim2
