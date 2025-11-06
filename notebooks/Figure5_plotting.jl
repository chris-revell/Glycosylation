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
using Statistics
using JLD2

@from "$(srcdir("Glycosylation.jl"))" using Glycosylation
@from "$(srcdir("Visualise.jl"))" using Visualise
@from "$(srcdir("UsefulFunctions.jl"))" using UsefulFunctions
@from "$(srcdir("MakeWeightMatrices.jl"))" using MakeWeightMatrices
@from "$(srcdir("DerivedParameters.jl"))" using DerivedParameters

#%%

subFolder = "Figure5"
folderName = "25-11-06-16-07-22_K₂=0.3_K₄=1.0_T̃ᵣ=3.85_differencing=centre_nSpatialDims=2_thicknessProfile=GRF_α_C=5.0_β=70.0_𝒟=204.0"
data1 = load(datadir("sims", subFolder, folderName, "solution.jld2"))
@unpack sol1, p1, sol2, p2, rawParams = data1
mat_h1 = reshape([p1.hᵥ[i,i] for i=1:prod(p1.dims)], p1.dims...)
@unpack thicknessProfile, differencing, solver, nOutputs, σGRF, λGRF, nSpatialDims, Ngrid, dims, h₀, 𝒜, Ω, N, k_Cd, k_Ca, k_Sd, k_Sa, k₁, k₂, k₃, k₄, 𝒞, 𝒮, ℰ, D_C, D_S, Tᵣstar, ϕ = rawParams
derivedParams = derivedParameters(Ω, 𝒜, N, k_Cd, k_Ca, k_Sd, k_Sa, k₁, k₂, k₃, k₄, 𝒞, 𝒮, ℰ, D_C, D_S, Tᵣstar; checks=false)
@unpack L₀, E₀, C_b, S_b, δ_C, δ_S, α_C, α_S, C₀, S₀, Tᵣ, T̃ᵣ, K₂, K₃, K₄, σ, ϵ, 𝒟, β, h_C, h_S, λ, ζ, γ, Δ, F = derivedParams

# uSlices = [selectdim(reshape(u, dims...), 3, dims[3]÷2) for u in sol1.u]
# uSlicesReshaped = [reshape(u, prod(dims[Not(3)])) for u in uSlices]
# concentrationSurfaceMovie(uSlicesReshaped, p1.dims[1:2]; subFolder=subFolder, folderName=folderName)
# concentrationHeatmapMovie(uSlicesReshaped, dims; subFolder=subFolder, folderName=folderName)
# M̃movie(sol1.u, p1; subFolder=subFolder, folderName=folderName)
# thicknessPlot(p1.hᵥ, p1.dims; subFolder=subFolder, folderName=folderName)


colorsUsed = [(:red), (:green), (:blue)]
νs = collect(range(0.0,1.0,p1.dims[1]))
xs = collect(range(0.0,sqrt(π),p1.dims[2]))

#%%

fig = Figure(size=(1000,800), fontsize=18, figure_padding=25)#, theme=textheme)
g1 = GridLayout(fig[1,1])
g2 = GridLayout(fig[2,1])


ax0 = Axis(g1[1, 1], aspect=DataAspect())
Label(g1[1,1,Bottom()], L"x")
Label(g1[1,1,Left()], L"y", rotation=π/2)
mat_h = reshape([p1.hᵥ[i,i] for i=1:prod(p1.dims)], p1.dims...)
maxdif = max(abs(minimum(mat_h)-1.0), abs(maximum(mat_h)-1.0))
clim = (1-maxdif, 1+maxdif)
heatmap!(ax0, collect(range(0.0, sqrt(π), p1.dims[1])), collect(range(0.0, sqrt(π), p1.dims[1])), mat_h[1,:,:], colorrange=clim)
ax0.xticks = (0.0:sqrt(π):sqrt(π), [L"0.0", L"\sqrt{\pi}"])
ax0.yticks = (0.0:sqrt(π):sqrt(π), [L"0.0", L"\sqrt{\pi}"])
text!(ax0, Point{2,Float64}(0.95*sqrt(π), 0.9*sqrt(π)), text="A", color=:white, align=(:right, :bottom), fontsize=24) 
indMax = findmax(mat_h1[1,:,p1.dims[3]÷2])[2]
indMin = findmin(mat_h1[1,:,p1.dims[3]÷2])[2]
peakxs = xs[[indMax, indMin]]
hlines!(ax0, sqrt(π)/2.0, color=(:white, 1.0), linewidth=2)
scatter!(ax0, peakxs, [sqrt(π)/2.0, sqrt(π)/2.0], marker=:star6, color=:white, markersize=20)
Colorbar(g1[1,2], limits=clim, label=L"h\left(\mathbf{x}_\perp\right)")


endPoint = length(sol2.u)-20
frameInds = collect(1:endPoint÷2-1:endPoint)

ax1 = CairoMakie.Axis(g1[1, 3])
uInternal = zeros(Float64, p1.dims[1:2]...)
for i in frameInds
    uInternal .= max.(uInternal, reshape(sol1.u[i], p1.dims...)[:,:,p1.dims[3]÷2])
end
Label(g1[1,3,Bottom()], L"\nu")
Label(g1[1,3,Left()], L"x", rotation=π/2)
ax1.yticks = (0.0:sqrt(π):sqrt(π), [L"0.0", L"\sqrt\pi"])
ax1.xticks = (0.0:1.0:1.0, [L"0.0", L"1.0"])
clim = (0.0, 30.0)
heatmap!(ax1, νs, xs, uInternal, colorrange=clim )
# scatter!(ax1, [0.0, 0.0], peakxs, marker=:star6, color=:white, markersize=20)
hlines!(ax1, peakxs, color=:white, linewidth=2)
text!(ax1, Point{2,Float64}(0.95,0.9*sqrt(π)), text="B", color=:white, align=(:right, :bottom), fontsize=24) 
Colorbar(g1[1,4], limits=clim, label=L"\tilde{C}\left(\nu, x, \sqrt{\pi}/2, \tilde{t}\right)")


ax2 = CairoMakie.Axis(g2[1, 1])
M̃local = zeros(Float64, p1.dims[1])
allLines = []
labels = []
times = []
for (c,i) in enumerate(frameInds)
    M̃local .= M̃(sol1.u[i], p1.W, p1.dims, p1.dν, p1.hᵥ)
    push!(allLines, lines!(ax2, collect(range(0.0, 1.0, p1.dims[1])), M̃local[:,1,1], linestyle=:solid, color=(colorsUsed[c], 0.5), linewidth=4))
    str = @sprintf("%.2f", sol1.t[i])
    push!(labels, L"\tilde{t}=%$(str),\ GRF")

    M̃local .= M̃(sol2.u[i], p2.W, p2.dims, p2.dν, p2.hᵥ)
    push!(allLines, lines!(ax2, collect(range(0.0, 1.0, p1.dims[1])), M̃local[:,1,1], linestyle=:dot, color=(colorsUsed[c], 1.0), linewidth=4))
    str = @sprintf("%.2f", sol1.t[i])
    push!(labels, L"\tilde{t}=%$(str),\ Uniform")
    push!(times, L"\tilde{t}=%$(str)")
end
# Label(g2[1,1,Bottom()], L"\nu")
Label(g2[1,1,Left()], L"\tilde{M}", rotation=π/2)
ax2.xticks = (0.0:1.0:1.0, [L"0.0", L"1.0"])
ax2.xticks = (0.0:0.5:1.0, [L"0.0", L"\phi", L"1.0"])
ax2.yticks = (0.0:50.0:50.0, [L"0.0", L"50.0"])    
ax2.xlabel = L"\nu"
vlines!(ax2, 0.5, color=(:black,0.5))
text!(ax2, Point{2,Float64}(0.95,0.9*50.0), text="C", color=:black, align=(:right, :bottom), fontsize=24) 
text!(ax2, Point{2,Float64}(0.05, 0.85*50.0), text = times[1], color=:red) 
text!(ax2, Point{2,Float64}(0.35, 0.45*50.0), text = times[2], color=:green) 
text!(ax2, Point{2,Float64}(0.7, 0.3*50.0), text = times[3], color=:blue) 
ax2.xgridvisible = false
ax2.ygridvisible = false
xlims!(ax2, (0.0, 1.0))
mlim = (0.0, 50.0)
ylims!(ax2, mlim)

colsize!(g1, 1, Aspect(1, 1.0))
colsize!(g1, 2, Aspect(1, 0.1))
colsize!(g1, 3, Aspect(1, 1.0))
colsize!(g1, 4, Aspect(1, 0.1))
colsize!(g2, 1, Aspect(1, 1.5))

resize_to_layout!(fig)
display(fig)
save(datadir("sims", subFolder, folderName, "Figure5.png"), fig)
save(datadir("sims", subFolder, folderName, "Figure5.pdf"), fig)

@show sol1.t[frameInds]
@show sol2.t[frameInds]

include(projectdir("notebooks", "paramsRaw.jl"))
derivedParams = derivedParameters(Ω, 𝒜, N, k_Cd, k_Ca, k_Sd, k_Sa, k₁, k₂, k₃, k₄, 𝒞, 𝒮, ℰ, D_C, D_S, Tᵣstar; checks=true)
@unpack L₀, E₀, C_b, S_b, δ_C, δ_S, α_C, α_S, C₀, S₀, Tᵣ, T̃ᵣ, K₂, K₃, K₄, σ, ϵ, 𝒟, β, h_C, h_S, λ, ζ, γ, Δ, F = derivedParams
ind50 = findfirst(x->M̃ϕ(x, p1.W, p1.dims, p1.dν, p1.hᵥ, 0.5) > 0.5*π, sol1.u)
Tᵣ₅₀Star = sol1.t[ind50]*(N^2)*(K₂+σ*K₃)/(k₁*E₀)    
𝒫sim1 = Mstarϕ(sol1.u[ind50], p1.W, p1.dims, p1.dν, p1.hᵥ, α_C, 𝒞, 0.5)/Tᵣ₅₀Star
@show 𝒫sim1

ind50 = findfirst(x->M̃ϕ(x, p2.W, p2.dims, p2.dν, p2.hᵥ, 0.5) > 0.5*π, sol2.u)
Tᵣ₅₀Star = sol2.t[ind50]*(N^2)*(K₂+σ*K₃)/(k₁*E₀)    
𝒫sim2 = Mstarϕ(sol2.u[ind50], p2.W, p2.dims, p2.dν, p2.hᵥ, α_C, 𝒞, 0.5)/Tᵣ₅₀Star
@show 𝒫sim2