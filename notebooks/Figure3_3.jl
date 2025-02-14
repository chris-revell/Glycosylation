
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
using LinearAlgebra

@from "$(srcdir("Glycosylation.jl"))" using Glycosylation
@from "$(srcdir("Visualise.jl"))" using Visualise
@from "$(srcdir("UsefulFunctions.jl"))" using UsefulFunctions
@from "$(srcdir("MakeWeightMatrices.jl"))" using MakeWeightMatrices
@from "$(srcdir("DerivedParameters.jl"))" using DerivedParameters
@from "$(srcdir("CisternaWidth.jl"))" using CisternaWidth

#%%

subFolder = "Figure3"
terminateAt = "nuWall"
# thicknessProfile = "Gaussian"
differencing = "centre"
solver = SSPRK432()
nOutputs = 100
# σGRF = 0.2
σGaussian = 0.20

nSpatialDims = 1
Ngrid = 401
dims = fill(Ngrid, nSpatialDims+1)

include(projectdir("notebooks", "paramsRaw.jl"))

#%%

rawParams1 = (
    thicknessProfile = "Gaussian",
    differencing = differencing,
    solver = solver,
    nOutputs = nOutputs,
    # σGRF = σGRF,
    nSpatialDims = nSpatialDims,
    Ngrid = Ngrid,
    dims = dims,
    h₀ = h₀,
    Ωperp = Ωperp,
    Ω = Ω,
    N = N,
    k_Cd = k_Cd,
    k_Ca = k_Ca,
    k_Sd = k_Sd,
    k_Sa = k_Sa,
    k₁ = k₁,
    k₂ = k₂,
    k₃ = k₃,
    k₄ = k₄,
    𝒞 = 𝒞,
    𝒮 = 𝒮,
    ℰ = ℰ,
    D_C = D_C,
    D_S = D_S,
    Tᵣstar = Tᵣstar,
    ϕ = ϕ
)

#%%

derivedParams = derivedParameters(Ω, Ωperp, N, k_Cd, k_Ca, k_Sd, k_Sa, k₁, k₂, k₃, k₄, 𝒞, 𝒮, ℰ, D_C, D_S, Tᵣstar; checks=true)
@unpack L₀, E₀, C_b, S_b, δ_C, δ_S, α_C, α_S, C₀, S₀, Tᵣ, T̃ᵣ, K₂, K₃, K₄, σ, ϵ, 𝒟, β, h_C, h_S, λ, ζ, γ, Δ, F = derivedParams

#%%

# Create directory for run data labelled with current time.
# paramsName = @savename nSpatialDims K₂ K₄ α_C β 𝒟 T̃ᵣ differencing
# folderName = "$(Dates.format(Dates.now(),"yy-mm-dd-HH-MM-SS"))_$(paramsName)"
# # Create frames subdirectory to store system state at each output time
# mkpath(datadir("sims",subFolder,folderName))

folderName = "25-02-07-14-18-04_K₂=0.3_K₄=1.0_T̃ᵣ=0.385_differencing=centre_nSpatialDims=1_α_C=5.0_β=70.0_𝒟=204.0"

#%%

# sol1, p1 = glycosylationAnyD(dims, K₂, K₄, T̃ᵣ, α_C, 𝒟, β, thickness="Gaussian", differencing=differencing, solver=solver, nOutputs=nOutputs, σGaussian=σGaussian, terminateAt=terminateAt)
# println("finished sim")
data1 = load(datadir("sims", subFolder, folderName, "solutionHVariation.jld2"))
@unpack sol1, p1 = data1

mat_h1 = reshape([p1.hᵥ[i,i] for i=1:prod(dims)], dims...)

#%%

# jldsave(datadir("sims",subFolder,folderName,"solutionHVariation.jld2"); sol1, p1, rawParams1)

#%%

# sol2, p2 = glycosylationAnyD(dims, K₂, K₄, T̃ᵣ, α_C, 𝒟, β, thickness="uniform", fDist="Gaussian", differencing=differencing, solver=solver, nOutputs=nOutputs, σGaussian=σGaussian, terminateAt=terminateAt)
# println("finished sim 2")
data2 = load(datadir("sims", subFolder, folderName, "solutionFVariation.jld2"))
@unpack sol2, p2 = data2

mat_h2 = reshape([p2.hᵥ[i,i] for i=1:prod(dims)], dims...)

#%%

rawParams2 = (
    thicknessProfile = "uniform",
    differencing = differencing,
    solver = solver,
    nOutputs = nOutputs,
    # σGRF = σGRF,
    nSpatialDims = nSpatialDims,
    Ngrid = Ngrid,
    dims = dims,
    h₀ = h₀,
    Ωperp = Ωperp,
    Ω = Ω,
    N = N,
    k_Cd = k_Cd,
    k_Ca = k_Ca,
    k_Sd = k_Sd,
    k_Sa = k_Sa,
    k₁ = k₁,
    k₂ = k₂,
    k₃ = k₃,
    k₄ = k₄,
    𝒞 = 𝒞,
    𝒮 = 𝒮,
    ℰ = ℰ,
    D_C = D_C,
    D_S = D_S,
    Tᵣstar = Tᵣstar,
    ϕ = ϕ
)

# jldsave(datadir("sims",subFolder,folderName,"solutionFVariation.jld2"); sol2, p2, rawParams2)

#%%

# frames = collect(4:48:100)
outLength = min(length(sol1.t), length(sol2.t))-5
frames = collect(1:outLength÷2-1:outLength)
# frames = collect(1:40:100)
fig = Figure(size=(1000,1000))

νs = collect(range(0.0, 1.0, dims[1]))
xMax = sqrt(π)
xs = collect(range(0.0, xMax, dims[2]))

letterArray = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)", "(i)", "(j)", "(k)", "(l)", "(m)", "(n)", "(o)", "(p)"]



g1 = GridLayout(fig[1,1])

ax1 = Axis(fig[1,2:4])
for x=1:3
    uInternal = reshape(sol1.u[frames[x]], dims...)
    heatmap!(ax1, νs.+(x-1)*1.0, xs, uInternal, colormap=:batlow)
end


display(fig)

Label(fig[2,x,Top()], popfirst!(letterArray))








axesVec = [Axis(fig[1,1], alignmode = Mixed(
    left = Makie.Protrusion(0),
    right = Makie.Protrusion(0),
    bottom = Makie.Protrusion(0),
    top = Makie.Protrusion(0),
))]
Label(fig[2,1,Top()], popfirst!(letterArray))
lines!(axesVec[1], mat_h1[1,:], xs)
xlims!(axesVec[1], (0.0, 1.2*maximum(mat_h1[1,:])))
ylims!(axesVec[1], (0.0, xMax))
axesVec[end].yticks = (0.0:sqrt(π)/2.0:sqrt(π), [L"0.0", L"\sqrt{\pi}/2", L"\sqrt{\pi}"])
for x=2:4
    uInternal = reshape(sol1.u[frames[x-1]], dims...)
    push!(axesVec, Axis(fig[1,x], alignmode = Mixed(
        left = Makie.Protrusion(0),
        right = Makie.Protrusion(0),
        bottom = Makie.Protrusion(0),
        top = Makie.Protrusion(0),
    )))
    # tString = @sprintf("%.3f", sol2.t[frames[x-1]])
    # Label(fig[1,x, Top()], L"\tilde{t} = %$tString")
    Label(fig[2,x,Top()], popfirst!(letterArray))
    heatmap!(axesVec[end], νs, xs, uInternal, colormap=:batlow)
end

for x=2:4  
    push!(axesVec, Axis(fig[3,x], alignmode = Mixed(
        left = Makie.Protrusion(0),
        right = Makie.Protrusion(0),
        bottom = Makie.Protrusion(0),
        top = Makie.Protrusion(0),
    )))
    Label(fig[4,x,Top()], popfirst!(letterArray))
    M = M̃(sol1.u[frames[x-1]], p1.W, dims, p1.dν, p1.hᵥ)[:,1]
    lines!(axesVec[end], νs, M)
end

push!(axesVec, Axis(fig[5,1], alignmode = Mixed(
    left = Makie.Protrusion(0),
    right = Makie.Protrusion(0),
    bottom = Makie.Protrusion(0),
    top = Makie.Protrusion(0),
)))
Label(fig[6,1,Top()], popfirst!(letterArray))
lines!(axesVec[end], p2.matFₑ, xs)
xlims!(axesVec[end], (0.0, 1.2*maximum(p2.matFₑ)))
ylims!(axesVec[end], (0.0, xMax))
axesVec[end].yticks = (0.0:sqrt(π)/2.0:sqrt(π), [L"0.0", L"\sqrt{\pi}/2", L"\sqrt{\pi}"])
for x=2:4
    uInternal = reshape(sol2.u[frames[x-1]], dims...)
    push!(axesVec, Axis(fig[5,x],alignmode = Mixed(
        left = Makie.Protrusion(0),
        right = Makie.Protrusion(0),
        bottom = Makie.Protrusion(0),
        top = Makie.Protrusion(0),
    )))
    Label(fig[6,x,Top()], popfirst!(letterArray))
    heatmap!(axesVec[end], νs, xs, uInternal, colormap=:batlow)
end

for x=2:4
    push!(axesVec, Axis(fig[7,x], alignmode = Mixed(
        left = Makie.Protrusion(0),
        right = Makie.Protrusion(0),
        bottom = Makie.Protrusion(0),
        top = Makie.Protrusion(0),
    )))
    Label(fig[8,x,Top()], popfirst!(letterArray))
    M = M̃(sol2.u[frames[x-1]], p2.W, dims, p2.dν, p2.hᵥ)[:,1]
    lines!(axesVec[end], νs, M)
end

# Colorbar(fig[3,5], limits=clim, label=L"\tilde{C}")

axesVec[1].xlabel = L"h"
axesVec[1].ylabel = L"x"
for ax in axesVec[2:4]
    # ax.xlabel = L"\nu"
    # ax.ylabel = L"x"
    ax.xticks = (0.0:0.5:1.0, [L"0.0", L"0.5", L"1.0"])
    ax.yticks = (0.0:sqrt(π)/2.0:sqrt(π), [L"0.0", L"\sqrt{\pi}/2", L"\sqrt{\pi}"])    
    ax.xticklabelsvisible = false
    ax.yticklabelsvisible = false
    xlims!(ax, (0.0,1.0))
    ylims!(ax, (0.0,xMax))
end
axesVec[5].ylabel = L"\tilde{M}"
for ax in axesVec[5:7]
    ax.xlabel = L"\nu"
    # ax.ylabel = L"\tilde{M}"
    xlims!(ax, (0.0,1.0))
    ylims!(ax, (0.0,40.0))
    ax.xticks = (0.0:0.5:1.0, [L"0.0", L"0.5", L"1.0"])
    ax.yticks = (0.0:20.0:40.0, [L"0.0", L"20.0", L"40.0"])
    # ax.xticklabelsvisible = false
    ax.yticklabelsvisible = false
end
axesVec[5].yticklabelsvisible = true

axesVec[8].xlabel = L"F_e"
axesVec[8].ylabel = L"x"
for ax in axesVec[9:11]
    # ax.xlabel = L"\nu"
    # ax.ylabel = L"x"
    ax.xticks = (0.0:0.5:1.0, [L"0.0", L"0.5", L"1.0"])
    ax.yticks = (0.0:sqrt(π)/2.0:sqrt(π), [L"0.0", L"\sqrt{\pi}/2", L"\sqrt{\pi}"])    
    ax.xticklabelsvisible = false
    ax.yticklabelsvisible = false
    xlims!(ax, (0.0,1.0))
    ylims!(ax, (0.0,xMax))
end
axesVec[12].ylabel = L"\tilde{M}"
for ax in axesVec[12:end]
    ax.xlabel = L"\nu"
    # ax.ylabel = L"\tilde{M}"
    xlims!(ax, (0.0,1.0))
    ylims!(ax, (0.0,40.0))
    ax.xticks = (0.0:0.5:1.0, [L"0.0", L"0.5", L"1.0"])
    ax.yticks = (0.0:20.0:40.0, [L"0.0", L"20.0", L"40.0"])
    # ax.xticklabelsvisible = false
    ax.yticklabelsvisible = false
end
axesVec[12].yticklabelsvisible = true

colsize!(fig.layout, 1, Aspect(1, 1.0))
colsize!(fig.layout, 2, Aspect(1, 1.0))
colsize!(fig.layout, 3, Aspect(1, 1.0))
colsize!(fig.layout, 4, Aspect(1, 1.0))
rowsize!(fig.layout, 2, Relative(0.002))
rowsize!(fig.layout, 4, Relative(0.002))
rowsize!(fig.layout, 6, Relative(0.002))
rowsize!(fig.layout, 8, Relative(0.002))

resize_to_layout!(fig)

display(fig)
save(datadir("sims",subFolder,folderName,"Figure3.png"), fig)
