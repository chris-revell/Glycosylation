
module Visualise

using LinearAlgebra
using SparseArrays
using CairoMakie
using DrWatson
using FromFile 

@from "$(srcdir("UsefulFunctions.jl"))" using UsefulFunctions

function concentrationSurfaceMovie(solu, ts, xs, νs, dims, Nghost, ghostVertexMaskVec; subFolder="", folderName="") 
    isdir(datadir("sims", subFolder, folderName)) ? nothing : mkdir(datadir("sims", subFolder, folderName))
    fig = Figure(size=(1000,1000))
    ax = Axis3(fig[1, 1], aspect=:equal, azimuth=-π/4)
    ax.xlabel = "ν"
    ax.ylabel = "x"
    ax.zlabel = "c"
    uInternal = Observable(zeros(dims))
    globalmin = minimum([minimum(u[ghostVertexMaskVec]) for u in solu])
    globalmax = maximum([maximum(u[ghostVertexMaskVec]) for u in solu])
    zlims!(ax, (globalmin, globalmax))
    clims = (globalmin,globalmax)
    surface!(ax, xs[Nghost+1:end-Nghost], νs[Nghost+1:end-Nghost], uInternal, colorrange=clims, colormap=:batlow)
    record(fig, datadir("sims",subFolder,folderName,"concentrationSurfaceMovie.mp4"), 1:length(ts); framerate=10) do i
        uInternal[] .= reshape(solu[i][ghostVertexMaskVec], dims)
        uInternal[] = uInternal[]
    end
    return nothing 
end

function spaceIntegralOver_ν_Movie(solu, ts, xs, νs, dims, Nghost, vertexWeightsMatrix, ghostVertexMaskVec; subFolder="", folderName="")
    isdir(datadir("sims", subFolder, folderName)) ? nothing : mkdir(datadir("sims", subFolder, folderName))
    # Find limits
    uInternal2D = reshape((vertexWeightsMatrix*solu[end])[ghostVertexMaskVec], dims)
    M = sum(uInternal2D, dims=2)[:,1]
    minima = Float64[]
    maxima = Float64[]
    for i=1:length(ts)
        uInternal2D .= reshape((vertexWeightsMatrix*solu[i])[ghostVertexMaskVec], dims)
        M .= sum(uInternal2D, dims=2)[:,1]
        push!(minima, minimum(M))
        push!(maxima, maximum(M))
    end
    globalmin = minimum(minima)
    globalmax = maximum(maxima)

    fig = Figure(size=(1000,1000))
    ax = CairoMakie.Axis(fig[1, 1], aspect=1)
    ax.xlabel = "ν"
    ax.ylabel = "M, ∱cdxdy"
    ax.title = "Integral of Cₛ over x against ν"
    M = Observable(zeros(dims[1]))
    lines!(ax, νs[1:Nghost:end-2*Nghost], M)
    ylims!(ax, (globalmin, globalmax))
    record(fig, datadir("sims",subFolder, folderName, "spaceIntegralOver_ν_Movie.mp4"), 1:length(ts); framerate=10) do i
        uInternal2D .= reshape((vertexWeightsMatrix*solu[i])[ghostVertexMaskVec], dims)
        M[] .= sum(uInternal2D, dims=2)[:,1]
        M[] = M[]
    end
    save(datadir("sims",subFolder,folderName,"finalSpaceIntegralOver_ν.png"), fig)
    return nothing
end

function productionHeatmap3D(ϕ, solu, ts, xs, νs, dims, ghostVertexMaskVec, vertexWeightsMatrix; subFolder="", folderName="")
    isdir(datadir("sims", subFolder, folderName)) ? nothing : mkdir(datadir("sims", subFolder, folderName))

    uInternal = reshape((vertexWeightsMatrix*solu[end])[ghostVertexMaskVec], dims)
    MsInternal = sum(uInternal[round(Int64, ϕ*dims[1]):end, :, :], dims=1)
    globalmax = maximum(MsInternal)

    fig = Figure(size=(1000,1000))
    ax = CairoMakie.Axis(fig[1, 1], aspect=1)
    ax.xlabel = "x"
    ax.ylabel = "y"
    ax.title = "Useful production P over x and y"
    M = Observable(zeros(dims[1:2]))
    
    heatmap!(ax, M, colorrange=(0.0, globalmax), colormap=:inferno)
    record(fig, datadir("sims",subFolder,folderName,"productionHeatmap.mp4"), 1:length(solu); framerate=10) do i
        uInternal .= reshape((vertexWeightsMatrix*solu[i])[ghostVertexMaskVec], dims)
        MsInternal .= sum(uInternal[round(Int64, ϕ*dims[1]):end, :, :], dims=1)
        M[] .= MsInternal[1,:,:]
        M[] = M[]
    end
    
    return nothing
end


export concentrationSurfaceMovie
export spaceIntegralOver_ν_Movie
export productionHeatmap3D

end
