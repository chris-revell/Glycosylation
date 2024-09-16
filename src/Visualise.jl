
module Visualise

using LinearAlgebra
using SparseArrays
using CairoMakie
using DrWatson
using FromFile 

@from "$(srcdir("UsefulFunctions.jl"))" using UsefulFunctions

# function concentrationSurfaceMovie(solu, ts, xs, νs, dimsReal, Nghost, ghostVertexMaskVec; subFolder="", folderName="") 
function concentrationSurfaceMovie(solu, ts, dimsReal, Nghost, ghostVertexMaskVec; subFolder="", folderName="") 
    isdir(datadir("sims", subFolder, folderName)) ? nothing : mkdir(datadir("sims", subFolder, folderName))
    fig = Figure(size=(1000,1000))
    ax = Axis3(fig[1, 1], aspect=:equal, azimuth=-π/4)
    ax.xlabel = "ν"
    ax.ylabel = "x"
    ax.zlabel = "c"
    uInternal = Observable(zeros(dimsReal...))
    globalmin = minimum([minimum(u[ghostVertexMaskVec]) for u in solu])
    globalmax = maximum([maximum(u[ghostVertexMaskVec]) for u in solu])
    zlims!(ax, (globalmin, globalmax))
    clims = (globalmin,globalmax)
    # surface!(ax, xs[Nghost+1:end-Nghost], νs[Nghost+1:end-Nghost], uInternal, colorrange=clims, colormap=:batlow)
    surface!(ax, uInternal, colorrange=clims, colormap=:batlow)
    record(fig, datadir("sims",subFolder,folderName,"concentrationSurfaceMovie.mp4"), 1:length(ts); framerate=10) do i
        uInternal[] .= reshape(solu[i][ghostVertexMaskVec], dimsReal...)
        uInternal[] = uInternal[]
    end
    return nothing 
end

function spaceIntegralOver_ν_Movie(solu, ts, xs, νs, dimsReal, Nghost, vertexWeightsMatrix, ghostVertexMaskVec; subFolder="", folderName="")
    isdir(datadir("sims", subFolder, folderName)) ? nothing : mkdir(datadir("sims", subFolder, folderName))
    dν = νs[3]-νs[2]
    # Find limits
    uReshaped = reshape((vertexWeightsMatrix*solu[end])[ghostVertexMaskVec], dimsReal...)
    M = sum(uReshaped, dims=(2,3))
    minima = Float64[]
    maxima = Float64[]
    for i=1:length(ts)
        uReshaped .= reshape((vertexWeightsMatrix*solu[i])[ghostVertexMaskVec], dimsReal...)
        M .= sum(uReshaped, dims=(2,3))./dν
        push!(minima, minimum(M))
        push!(maxima, maximum(M))
    end
    globalmin = minimum(minima)
    globalmax = maximum(maxima)

    fig = Figure(size=(1000,1000))
    ax = CairoMakie.Axis(fig[1, 1], aspect=1)
    ax.xlabel = "ν"
    ax.ylabel = "M, ∱cdxdy"
    ax.title = "Integral of Cₛ over space against ν"
    M = Observable(zeros(dimsReal[1]))
    lines!(ax, νs[2:end-1], M)
    ylims!(ax, (globalmin, globalmax))
    record(fig, datadir("sims",subFolder, folderName, "spaceIntegralOver_ν_Movie.mp4"), 1:length(ts); framerate=10) do i
        uReshaped .= reshape((vertexWeightsMatrix*solu[i])[ghostVertexMaskVec], dimsReal...)
        M[] .= dropdims((sum(uReshaped, dims=(2,3))), dims=Tuple(collect(2:ndims(uReshaped))))./dν
        M[] = M[]
    end
    save(datadir("sims",subFolder,folderName,"finalSpaceIntegralOver_ν.png"), fig)
    return nothing
end

function productionHeatmap3D(ϕ, solu, ts, xs, νs, dimsReal, ghostVertexMaskVec, vertexWeightsMatrix; subFolder="", folderName="")
    isdir(datadir("sims", subFolder, folderName)) ? nothing : mkdir(datadir("sims", subFolder, folderName))

    uInternal = reshape((vertexWeightsMatrix*solu[end])[ghostVertexMaskVec], dimsReal...)
    MsInternal = sum(uInternal[round(Int64, ϕ*dimsReal[1]):end, :, :], dims=1)
    globalmax = maximum(MsInternal)

    fig = Figure(size=(1000,1000))
    ax = CairoMakie.Axis(fig[1, 1], aspect=1)
    ax.xlabel = "x"
    ax.ylabel = "y"
    ax.title = "Useful production P over x and y"
    M = Observable(zeros(dimsReal[1:2]))
    
    heatmap!(ax, M, colorrange=(0.0, globalmax), colormap=:inferno)
    record(fig, datadir("sims",subFolder,folderName,"productionHeatmap.mp4"), 1:length(solu); framerate=10) do i
        uInternal .= reshape((vertexWeightsMatrix*solu[i])[ghostVertexMaskVec], dimsReal)
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
