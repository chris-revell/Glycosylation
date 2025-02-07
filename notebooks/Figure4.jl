using MathTeXEngine # required for texfont

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

#%%

subFolder = "Figure4"
terminateAt = "nuWall"
thicknessProfile = "GRF"
differencing = "centre"
solver = SSPRK432()
nOutputs = 100
σGRF = 0.3
λGRF = 0.2

nSpatialDims = 2
Ngrid = 101
dims = fill(Ngrid, nSpatialDims+1)

include(projectdir("notebooks", "paramsRaw.jl"))

#%%

derivedParams = derivedParameters(Ω, Ωperp, N, k_Cd, k_Ca, k_Sd, k_Sa, k₁, k₂, k₃, k₄, 𝒞, 𝒮, ℰ, D_C, D_S, Tᵣstar; checks=true)
@unpack L₀, E₀, C_b, S_b, δ_C, δ_S, α_C, α_S, C₀, S₀, Tᵣ, T̃ᵣ, K₂, K₃, K₄, σ, ϵ, 𝒟, β, h_C, h_S, λ, ζ, γ, Δ, F = derivedParams

#%%

# Create directory for run data labelled with current time.
paramsName = @savename nSpatialDims K₂ K₄ α_C β 𝒟 T̃ᵣ thicknessProfile differencing
folderName = "$(Dates.format(Dates.now(),"yy-mm-dd-HH-MM-SS"))_$(paramsName)"
# Create frames subdirectory to store system state at each output time
mkpath(datadir("sims",subFolder,folderName))

#%%

sol1, p1 = glycosylationAnyD(dims, K₂, K₄, T̃ᵣ, α_C, 𝒟, β, thickness=thicknessProfile, differencing=differencing, solver=solver, nOutputs=nOutputs, σGRF=σGRF, λGRF=λGRF, terminateAt=terminateAt)
println("finished sim 1")

sol2, p2 = glycosylationAnyD(fill(Ngrid, 2), K₂, K₄, T̃ᵣ, α_C, 𝒟, β, thickness="uniform", differencing=differencing, solver=solver, nOutputs=nOutputs, terminateAt=terminateAt)
println("finished sim 1")

#%%

rawParams = (
    thicknessProfile = thicknessProfile,
    differencing = differencing,
    solver = solver,
    nOutputs = nOutputs,
    σGRF = σGRF,
    λGRF = λGRF,
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

jldsave(datadir("sims",subFolder,folderName,"solution.jld2"); sol1, p1, sol2, p2, rawParams)

textheme = Theme(fonts=(; regular=texfont(:text),
                        bold=texfont(:bold),
                        italic=texfont(:italic),
                        bold_italic=texfont(:bolditalic)),
                        fontsize=32)

#%%

# fig = Figure(size=(1000,1000), theme=textheme)
# ax = CairoMakie.Axis(fig[1, 1])

# allLines = []
# labels = []
colorsUsed = [(:red), (:green), (:blue)]
# for (i,u) in enumerate(sol1.u[1:length(sol1.u)÷2-1:end])
#     M̃local = M̃(u, p.W, p.dims, p.dν, p.hᵥ)
#     # uInternal = reshape(p.W*p.hᵥ*u, p.dims...)
#     # M̃local = sum(uInternal, dims=(2:length(dims)))./dν
#     push!(allLines, lines!(ax, collect(range(0.0, 1.0, p.dims[1])), M̃local[:,1,1], linestyle=:solid, color=(colorsUsed[i], 0.5), linewidth=4))
#     str = @sprintf("%.2f", sol1.t[i])
#     push!(labels, L"\tilde{t}=%$(str)")
# end

# axislegend(ax, allLines, labels)
# ax.xlabel = L"\nu"
# ax.ylabel = L"\tilde{M}"
# ylims!(ax, (0.0, maximum(sol1.u[1])))
# xlims!(ax, (0.0, 1.0))

# save(datadir("sims", subFolder, folderName, "M̃.png"), fig)
# # display(fig)


# #%%

# fig = Figure(size=(1000,1000), theme=textheme)
# allAxes = []
# # allLines = []
# # labels = []
# # colorsUsed = [(:red), (:green), (:blue)]
# for (i,u) in enumerate(sol1.u[1:length(sol1.u)÷2-1:end])
#     push!(allAxes, CairoMakie.Axis(fig[i, 1], aspect=1))
#     uInternal = reshape(u, p.dims...)[:,:,p.dims[3]÷2]
#     allAxes[end].xlabel = L"\nu"
#     allAxes[end].ylabel = L"x"
#     # globalmin = minimum([minimum(u) for u in solu])
#     # globalmax = maximum([maximum(u) for u in solu])
#     # clims = (globalmin,globalmax)
#     νs = collect(range(0,1,p.dims[1]))
#     xs = collect(range(0,sqrt(π),p.dims[2]))
#     heatmap!(allAxes[end], νs, xs, uInternal)#, colorrange=clims, colormap=:batlow)

#     push!(allAxes, CairoMakie.Axis(fig[i, 2]))    
#     M̃local = M̃(u, p.W, p.dims, p.dν, p.hᵥ)
#     lines!(allAxes[end], collect(range(0.0, 1.0, p.dims[1])), M̃local[:,1,1], linestyle=:solid, color=(colorsUsed[i], 0.5), linewidth=4)
#     # str = @sprintf("%.2f", sol1.t[i])
#     # push!(labels, L"\tilde{t}=%$(str)")
#     allAxes[end].xlabel = L"\nu"
#     allAxes[end].ylabel = L"\tilde{M}"
# end

# display(fig)
# save(datadir("sims", subFolder, folderName, "Figure6.png"), fig)

# thicknessPlot(p.hᵥ, p.dims; subFolder=subFolder, folderName=folderName) 


#%%

fig = Figure(size=(1500,500), figure_padding=30, theme=textheme)

mat_h = reshape([p1.hᵥ[i,i] for i=1:prod(p1.dims)], p1.dims...)

ax0 = Axis(fig[1, 1], aspect=DataAspect())
ax0.xlabel = L"x"
ax0.ylabel = L"y"
maxdif = max(abs(minimum(mat_h)-1.0), abs(maximum(mat_h)-1.0))
clim = (1-maxdif, 1+maxdif)
heatmap!(ax0, collect(range(0.0, sqrt(π), p1.dims[1])), collect(range(0.0, sqrt(π), p1.dims[1])), mat_h[1,:,:], colorrange=clim)
ax0.xticks = (0.0:sqrt(π):sqrt(π), [L"0.0", L"\sqrt{\pi}"])
ax0.yticks = (0.0:sqrt(π):sqrt(π), [L"0.0", L"\sqrt{\pi}"])

indMax = findmax(mat_h[1,:,p1.dims[3]÷2])[2]
indMin = findmin(mat_h[1,:,p1.dims[3]÷2])[2]
peakxs = xs[[indMax, indMin]]
hlines!(ax0, sqrt(π)/2.0, color=(:black, 0.5), linewidth=4)
scatter!(ax0, peakxs, [sqrt(π)/2.0, sqrt(π)/2.0], marker=:star6, color=:orange, markersize=20)

Colorbar(fig[1,2], limits=clim, label=L"h(x)")



ax1 = CairoMakie.Axis(fig[1, 3])
uInternal = zeros(Float64, p1.dims[1:2]...)
for (i, u) in enumerate(sol1.u[1:length(sol1.u)÷2-1:end])
    uInternal .= max.(uInternal, reshape(u, p1.dims...)[:,:,p1.dims[3]÷2])
end
ax1.xlabel = L"\nu"
ax1.ylabel = L"x"
# globalmin = minimum([minimum(u) for u in solu])
# globalmax = maximum([maximum(u) for u in solu])
# clims = (globalmin,globalmax)
νs = collect(range(0.0,1.0,p1.dims[1]))
xs = collect(range(0.0,sqrt(π),p1.dims[2]))
ax1.yticks = (0.0:sqrt(π):sqrt(π), [L"0.0", L"\sqrt\pi"])
clim = (0.0, 30.0)
heatmap!(ax1, νs, xs, uInternal, colorrange=clim )#, colorrange=clims, colormap=:batlow)

scatter!(ax1, [0.0, 0.0], peakxs, marker=:star6, color=:orange, markersize=20)
hlines!(ax1, peakxs, color=:orange, linewidth=4)

Colorbar(fig[1,4], limits=clim, label=L"\tilde{C}(x)")

ax2 = CairoMakie.Axis(fig[1, 5])
M̃local = zeros(Float64, p1.dims[1])
allLines = []
labels = []
for (c,i) in enumerate([1, min(length(sol1.t), length(sol2.t))÷2, min(length(sol1.t), length(sol2.t))])
    M̃local .= M̃(sol1.u[i], p1.W, p1.dims, p1.dν, p1.hᵥ)
    push!(allLines, lines!(ax2, collect(range(0.0, 1.0, p1.dims[1])), M̃local[:,1,1], linestyle=:solid, color=(colorsUsed[c], 0.5), linewidth=4))
    str = @sprintf("%.2f", sol1.t[i])
    push!(labels, L"\tilde{t}=%$(str), GRF")

    M̃local .= M̃(sol2.u[i], p2.W, p2.dims, p2.dν, p2.hᵥ)
    push!(allLines, lines!(ax2, collect(range(0.0, 1.0, p1.dims[1])), M̃local[:,1,1], linestyle=:dot, color=(colorsUsed[c], 0.5), linewidth=4))
    str = @sprintf("%.2f", sol1.t[i])
    push!(labels, L"\tilde{t}=%$(str), Uniform")
end
ax2.xlabel = L"\nu"
ax2.ylabel = L"\tilde{M}"
# axislegend(ax2, allLines, labels)
Legend(fig[1,6], allLines, labels)
xlims!(ax2, (0.0, 1.0))
mlim = (0.0, 50.0)
ylims!(ax2, mlim)#maximum(M̃(sol1.u[1], p1.W, p1.dims, p1.dν, p1.hᵥ))))

tight_ticklabel_spacing!(ax0)


colsize!(fig.layout, 1, Aspect(1, 1.0))
colsize!(fig.layout, 3, Aspect(1, 1.0))
colsize!(fig.layout, 5, Aspect(1, 1.0))

Label(fig[2,1,Top()], L"(a)")
Label(fig[2,3,Top()], L"(b)")
Label(fig[2,5,Top()], L"(c)")
rowsize!(fig.layout, 2, Relative(0.02))

resize_to_layout!(fig)
display(fig)
save(datadir("sims", subFolder, folderName, "Figure4.png"), fig)


