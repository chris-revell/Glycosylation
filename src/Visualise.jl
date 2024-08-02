
module Visualise

using LinearAlgebra
using SparseArrays
using CairoMakie
using DrWatson

function concentrationSurfaceMovie(solu, ts, xs, νs, dims, Nghost, ghostVertexMask; subFolder="", folderName="") 
    isdir(datadir("sims", subFolder, folderName)) ? nothing : mkdir(datadir("sims", subFolder, folderName))
    fig = Figure(size=(1000,1000))
    ax = Axis3(fig[1, 1], aspect=:equal, azimuth=-π/4)
    ax.xlabel = "ν"
    ax.ylabel = "x"
    ax.zlabel = "c"
    uInternal = Observable(zeros(dims))
    globalmin = minimum([minimum(u[ghostVertexMask]) for u in solu])
    globalmax = maximum([maximum(u[ghostVertexMask]) for u in solu])
    zlims!(ax, (globalmin, globalmax))
    clims = (globalmin,globalmax)
    surface!(ax, xs[Nghost+1:end-Nghost], νs[Nghost+1:end-Nghost], uInternal, colorrange=clims, colormap=:batlow)
    record(fig, datadir("sims",subFolder,folderName,"concentrationSurfaceMovie.mp4"), 1:length(ts); framerate=10) do i
        uInternal[] .= reshape(solu[i][ghostVertexMask], dims)
        uInternal[] = uInternal[]
    end
    return nothing 
end

function spaceIntegralOver_ν_Movie(solu, ts, xs, νs, dims, Nghost, vertexWeightsMatrix, ghostVertexMask; subFolder="", folderName="")
    isdir(datadir("sims", subFolder, folderName)) ? nothing : mkdir(datadir("sims", subFolder, folderName))
    # Find limits
    uInternal2D = reshape((vertexWeightsMatrix*solu[end])[ghostVertexMask], dims)
    M = sum(uInternal2D, dims=2)[:,1]
    minima = Float64[]
    maxima = Float64[]
    for i=1:length(ts)
        uInternal2D .= reshape((vertexWeightsMatrix*solu[i])[ghostVertexMask], dims)
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
        uInternal2D .= reshape((vertexWeightsMatrix*solu[i])[ghostVertexMask], dims)
        M[] .= sum(uInternal2D, dims=2)[:,1]
        M[] = M[]
    end
    save(datadir("sims",subFolder,folderName,"finalSpaceIntegralOver_ν.png"), fig)
    return nothing
end


export concentrationSurfaceMovie
export spaceIntegralOver_ν_Movie

end


# fig = Figure(size=(1000,1000))
# ax = CairoMakie.Axis(fig[1, 1], aspect=1)
# ax.xlabel = "x"
# ax.ylabel = "y"
# uInternal3D = reshape(sol.u[1][ghostVertexMask], (Nx, Ny, Nν))
# uInternal = Observable(zeros(Nx,Ny))
# globalmin = minimum([minimum(u) for u in sol.u])
# globalmax = maximum([maximum(u) for u in sol.u])
# clims = (globalmin,globalmax)
# clims = (minimum(uInternal3D), maximum(uInternal3D))
# heatmap!(ax, xs[Nghost+1:end-Nghost], ys[Nghost+1:end-Nghost], uInternal, colorrange=clims, colormap=:batlow)
# record(fig, datadir("sims",folderName, "LargeNuTimeScan.mp4"), 1:length(ts); framerate=10) do i
#     uInternal3D .= reshape(sol.u[i][ghostVertexMask], (Nx, Ny, Nν))
#     uInternal[] .= uInternal3D[:,:,end]
#     uInternal[] = uInternal[]
#     ax.title = "xy profile of ν=1.0 at t=$(ts[i])"
# end

# #%%

# fig = Figure(size=(1000,1000))
# ax = CairoMakie.Axis(fig[1, 1], aspect=1)
# ax.xlabel = "ν"
# ax.ylabel = "c"
# uInternal3D = reshape(sol.u[1][ghostVertexMask], (Nx, Ny, Nν))
# uInternal = Observable(zeros(Nν))
# globalmin = minimum([minimum(u) for u in sol.u])
# globalmax = maximum([maximum(u) for u in sol.u])
# ylims!(ax, (globalmin, globalmax))
# lines!(ax, νs[Nghost+1:end-Nghost], uInternal)
# ax.title = "c against ν at x=0.5, y=0.5 at t=0.0"
# record(fig, datadir("sims",folderName, "NuProfileAtxyOverTime.mp4"), 1:length(ts); framerate=10) do i
#     uInternal3D .= reshape(sol.u[i][ghostVertexMask], (Nx, Ny, Nν))
#     uInternal[] .= uInternal3D[Nx÷2,Nx÷2,:]
#     uInternal[] = uInternal[]
#     ax.title = "c against ν at x=0.5, y=0.5 at t=$(ts[i])"
# end

# #%%

# fig = Figure(size=(1000,1000))
# ax = CairoMakie.Axis(fig[1, 1], aspect=1)
# ax.xlabel = "x"
# ax.ylabel = "ν"
# globalmin = minimum([minimum(u) for u in sol.u])
# globalmax = maximum([maximum(u) for u in sol.u])
# uInternal = Observable(zeros(Nx,Nν))
# clims = (globalmin,globalmax)
# ax.title = "c against x and ν at y=0.5 at t=0.0"
# heatmap!(ax, xs[Nghost+1:end-Nghost], νs[Nghost+1:end-Nghost], uInternal, colorrange=clims, colormap=:batlow)
# record(fig, datadir("sims",folderName, "xνOverTimeAty.mp4"), 1:length(ts); framerate=10) do i
#     uInternal[] .= reshape(sol.u[i][ghostVertexMask], (Nx, Ny, Nν))[:,Nyplus÷2,:]
#     uInternal[] = uInternal[]
#     ax.title = "c against x and ν at y=0.5 at t=$(ts[i])"
# end