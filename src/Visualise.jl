
module Visualise

using LinearAlgebra
using SparseArrays
using CairoMakie
using DrWatson
using FromFile 

@from "$(srcdir("UsefulFunctions.jl"))" using UsefulFunctions

function thicknessPlot(hᵥ, dims; subFolder="", folderName="") 
    isdir(datadir("sims", subFolder, folderName)) ? nothing : mkpath(datadir("sims", subFolder, folderName))
   
    mat_h = reshape([hᵥ[i,i] for i=1:prod(dims)], dims...)

    if length(dims) == 2
        fig = Figure()#size=(600,500))
        ax = Axis(fig[1, 1])    
        ax.xlabel = L"x"
        ax.ylabel = L"h"
        lines!(ax, mat_h[1,:])
        save(datadir("sims",subFolder,folderName,"thicknessPlot.png"), fig)
    else
        fig = Figure()#size=(600,500))
        ax = Axis(fig[1, 1], aspect=DataAspect())
        ax.xlabel = L"x"
        ax.ylabel = L"y"
        heatmap!(ax, mat_h[1,:,:])
        maxdif = max(abs(minimum(mat_h)-1.0), abs(maximum(mat_h)-1.0))
        Colorbar(fig[1,2], limits=(1-maxdif, 1+maxdif), label=L"h(x)")
        save(datadir("sims",subFolder,folderName,"thicknessPlot.png"), fig)
    end
    return nothing 
end

function concentrationSurfaceMovie(solu, dims; subFolder="", folderName="") 
    isdir(datadir("sims", subFolder, folderName)) ? nothing : mkpath(datadir("sims", subFolder, folderName))
    fig = Figure(size=(600,450))
    ax = Axis3(fig[1, 1], aspect=:equal, azimuth=-π/4)
    ax.xlabel = L"\nu"
    ax.ylabel = L"x"
    ax.zlabel = L"\tilde{C}(x_\perp, \nu)"
    uInternal = Observable(zeros(dims...))
    globalmin = minimum([minimum(u) for u in solu])
    globalmax = maximum([maximum(u) for u in solu])
    zlims!(ax, (globalmin, globalmax))
    clims = (globalmin,globalmax)
    # hidedecorations!(ax)
    νs = collect(range(0,1,dims[1]))
    xs = collect(range(0,sqrt(π),dims[2]))
    surface!(ax, νs, xs, uInternal, colorrange=clims, colormap=:batlow)
    record(fig, datadir("sims",subFolder,folderName,"concentrationSurfaceMovie.mp4"), 1:length(solu); framerate=10) do i
        uInternal[] .= reshape(solu[i], dims...)
        uInternal[] = uInternal[]
    end
    return nothing 
end

function concentrationHeatmapMovie(solu, dims; subFolder="", folderName="") 
    isdir(datadir("sims", subFolder, folderName)) ? nothing : mkpath(datadir("sims", subFolder, folderName))
    fig = Figure(size=(500,500))
    # ax = Axis(fig[1, 1], aspect=:equal)
    ax = Axis(fig[1, 1])
    ax.xlabel = L"\nu"
    ax.ylabel = L"x"
    uInternal = Observable(zeros(dims...))
    globalmin = minimum([minimum(u) for u in solu])
    globalmax = maximum([maximum(u) for u in solu])
    clims = (globalmin,globalmax)
    νs = collect(range(0,1,dims[1]))
    xs = collect(range(0,sqrt(π),dims[2]))
    heatmap!(ax, νs, xs, uInternal, colorrange=clims, colormap=:batlow)
    record(fig, datadir("sims",subFolder,folderName,"concentrationHeatmapMovie.mp4"), 1:length(solu); framerate=10) do i
        uInternal[] .= reshape(solu[i], dims...)
        uInternal[] = uInternal[]
    end
    return nothing 
end

function M̃movie(solu, p; subFolder="", folderName="")
    isdir(datadir("sims", subFolder, folderName)) ? nothing : mkpath(datadir("sims", subFolder, folderName))
    @unpack dims, dν, W, hᵥ = p
    # Find limits
    M = M̃(solu[end], W, dims, dν, hᵥ)
    minima = Float64[]
    maxima = Float64[]
    for u in solu        
        M .= M̃(u, W, dims, dν, hᵥ)
        push!(minima, minimum(M))
        push!(maxima, maximum(M))
    end
    globalmin = minimum(minima)
    globalmax = maximum(maxima)

    fig = Figure()#size=(500,500))
    ax = CairoMakie.Axis(fig[1, 1], aspect=1)
    ax.xlabel = L"\nu"
    ax.ylabel = L"\tilde{M}"
    # ax.title = "Integral of Cₛ over space against ν"
    M = Observable(zeros(dims[1]))
    lines!(ax, collect(range(0.0,1.0,dims[1])), M, linewidth=4)
    ylims!(ax, (globalmin, globalmax))
    record(fig, datadir("sims",subFolder, folderName, "Mtildemovie.mp4"), 1:length(solu); framerate=10) do i
        M[] .= M̃(solu[i], W, dims, dν, hᵥ)
        M[] = M[]
    end
    save(datadir("sims",subFolder,folderName,"Mtildefinal.png"), fig)
    return nothing
end

# function productionHeatmap3D(ϕ, solu, ts, xs, νs, dims, W; subFolder="", folderName="")
#     isdir(datadir("sims", subFolder, folderName)) ? nothing : mkpath(datadir("sims", subFolder, folderName))

#     uInternal = reshape((W*solu[end]), dims...)
#     MsInternal = sum(uInternal[round(Int64, ϕ*dims[1]):end, :, :], dims=1)
#     globalmax = maximum(MsInternal)

#     fig = Figure(size=(1000,1000))
#     ax = CairoMakie.Axis(fig[1, 1], aspect=1)
#     ax.xlabel = "x"
#     ax.ylabel = "y"
#     ax.title = "Useful production P over x and y"
#     M = Observable(zeros(dims[2:3]))
#     heatmap!(ax, M, colorrange=(0.0, globalmax), colormap=:inferno)
#     record(fig, datadir("sims",subFolder,folderName,"productionHeatmap.mp4"), 1:length(solu); framerate=10) do i
#         uInternal .= reshape((W*solu[i]), dims)
#         MsInternal .= sum(uInternal[round(Int64, ϕ*dims[1]):end, :, :], dims=1)
#         M[] .= MsInternal[1,:,:]
#         M[] = M[]
#     end
    
#     return nothing
# end


export concentrationSurfaceMovie
export concentrationHeatmapMovie
export M̃movie
# export productionHeatmap3D
export thicknessPlot

end
