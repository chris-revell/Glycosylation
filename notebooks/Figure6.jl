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
using MathTeXEngine # required for texfont

@from "$(srcdir("UsefulFunctions.jl"))" using UsefulFunctions
@from "$(srcdir("Visualise.jl"))" using Visualise

subFolder = "2spatialD"
folderName = "25-01-30-12-49-50_K₂=0.16_K₄=1.0_T̃ᵣ=0.431_differencing=centre_nSpatialDims=2_thicknessProfile=GRF_α_C=3.0_β=84.0_𝒟=1090.0"
# Create frames subdirectory to store system state at each output time
data = load(datadir("sims", subFolder, folderName, "solution.jld2"))
@unpack sol, p = data

textheme = Theme(fonts=(; regular=texfont(:text),
                        bold=texfont(:bold),
                        italic=texfont(:italic),
                        bold_italic=texfont(:bolditalic)),
                        fontsize=32)

#%%

fig = Figure(size=(1000,1000), theme=textheme)
ax = CairoMakie.Axis(fig[1, 1])

allLines = []
labels = []
colorsUsed = [(:red), (:green), (:blue)]
for (i,u) in enumerate(sol.u[1:length(sol.u)÷2-1:end])
    M̃local = M̃(u, p.W, p.dims, p.dν, p.hᵥ)
    # uInternal = reshape(p.W*p.hᵥ*u, p.dims...)
    # M̃local = sum(uInternal, dims=(2:length(dims)))./dν
    push!(allLines, lines!(ax, collect(range(0.0, 1.0, p.dims[1])), M̃local[:,1,1], linestyle=:solid, color=(colorsUsed[i], 0.5), linewidth=4))
    str = @sprintf("%.2f", sol.t[i])
    push!(labels, L"\tilde{t}=%$(str)")
end

axislegend(ax, allLines, labels)
ax.xlabel = L"\nu"
ax.ylabel = L"\tilde{M}"
ylims!(ax, (0.0, maximum(sol.u[1])))
xlims!(ax, (0.0, 1.0))

save(datadir("sims", subFolder, folderName, "M̃.png"), fig)
# display(fig)


#%%

fig = Figure(size=(1000,1000), theme=textheme)
allAxes = []
# allLines = []
# labels = []
# colorsUsed = [(:red), (:green), (:blue)]
for (i,u) in enumerate(sol.u[1:length(sol.u)÷2-1:end])
    push!(allAxes, CairoMakie.Axis(fig[i, 1], aspect=1))
    uInternal = reshape(u, p.dims...)[:,:,p.dims[3]÷2]
    allAxes[end].xlabel = L"\nu"
    allAxes[end].ylabel = L"x"
    # globalmin = minimum([minimum(u) for u in solu])
    # globalmax = maximum([maximum(u) for u in solu])
    # clims = (globalmin,globalmax)
    νs = collect(range(0,1,p.dims[1]))
    xs = collect(range(0,sqrt(π),p.dims[2]))
    heatmap!(allAxes[end], νs, xs, uInternal)#, colorrange=clims, colormap=:batlow)

    push!(allAxes, CairoMakie.Axis(fig[i, 2]))    
    M̃local = M̃(u, p.W, p.dims, p.dν, p.hᵥ)
    lines!(allAxes[end], collect(range(0.0, 1.0, p.dims[1])), M̃local[:,1,1], linestyle=:solid, color=(colorsUsed[i], 0.5), linewidth=4)
    # str = @sprintf("%.2f", sol.t[i])
    # push!(labels, L"\tilde{t}=%$(str)")
    allAxes[end].xlabel = L"\nu"
    allAxes[end].ylabel = L"\tilde{M}"
end

display(fig)
save(datadir("sims", subFolder, folderName, "Figure6.png"), fig)

thicknessPlot(p.hᵥ, p.dims; subFolder=subFolder, folderName=folderName) 


#%%

fig = Figure(size=(1500,500), figure_padding=30, theme=textheme)

mat_h = reshape([p.hᵥ[i,i] for i=1:prod(p.dims)], p.dims...)

ax0 = Axis(fig[1, 1], aspect=DataAspect())
ax0.xlabel = L"x"
ax0.ylabel = L"y"
heatmap!(ax0, collect(range(0.0, sqrt(π), p.dims[1])), collect(range(0.0, sqrt(π), p.dims[1])), mat_h[1,:,:])
ax0.xticks = (0.0:sqrt(π):sqrt(π), [L"0.0", L"\sqrt{\pi}"])
ax0.yticks = (0.0:sqrt(π):sqrt(π), [L"0.0", L"\sqrt{\pi}"])
maxdif = max(abs(minimum(mat_h)-1.0), abs(maximum(mat_h)-1.0))
Colorbar(fig[1,2], limits=(1-maxdif, 1+maxdif), label=L"h(x)")

ax1 = CairoMakie.Axis(fig[1, 3])
uInternal = zeros(Float64, p.dims[1:2]...)
for (i, u) in enumerate(sol.u[1:length(sol.u)÷2-1:end])
    uInternal .= max.(uInternal, reshape(u, p.dims...)[:,:,p.dims[3]÷2])
end
ax1.xlabel = L"\nu"
ax1.ylabel = L"x"
# globalmin = minimum([minimum(u) for u in solu])
# globalmax = maximum([maximum(u) for u in solu])
# clims = (globalmin,globalmax)
νs = collect(range(0.0,1.0,p.dims[1]))
xs = collect(range(0.0,sqrt(π),p.dims[2]))
ax1.yticks = (0.0:sqrt(π):sqrt(π), [L"0.0", L"\sqrt\pi"])
clim = (0.0, 30.0)
heatmap!(ax1, νs, xs, uInternal, colorrange=clim )#, colorrange=clims, colormap=:batlow)
Colorbar(fig[1,4], limits=clim, label=L"\tilde{C}(x)")

ax2 = CairoMakie.Axis(fig[1, 5])
M̃local = zeros(Float64, p.dims[1])
allLines = []
labels = []
for (i, u) in enumerate(sol.u[1:length(sol.u)÷2-1:end])
    M̃local = M̃(u, p.W, p.dims, p.dν, p.hᵥ)
    push!(allLines, lines!(ax2, collect(range(0.0, 1.0, p.dims[1])), M̃local[:,1,1], linestyle=:solid, color=(colorsUsed[i], 0.5), linewidth=4))
    str = @sprintf("%.2f", sol.t[i])
    push!(labels, L"\tilde{t}=%$(str)")
end
ax2.xlabel = L"\nu"
ax2.ylabel = L"\tilde{M}"
axislegend(ax2, allLines, labels)
xlims!(ax2, (0.0, 1.0))
mlim = (0.0, 50.0)
ylims!(ax2, mlim)#maximum(M̃(sol.u[1], p.W, p.dims, p.dν, p.hᵥ))))

tight_ticklabel_spacing!(ax0)


colsize!(fig.layout, 1, Aspect(1, 1.0))
colsize!(fig.layout, 3, Aspect(1, 1.0))
colsize!(fig.layout, 5, Aspect(1, 1.0))

# Label(fig[2,1,Top()], L"(a)")
# Label(fig[2,3,Top()], L"(b)")
# Label(fig[2,5,Top()], L"(c)")

# rowsize!(fig.layout, 2, Relative(0.2))


resize_to_layout!(fig)
display(fig)
save(datadir("sims", subFolder, folderName, "Figure6.png"), fig)


