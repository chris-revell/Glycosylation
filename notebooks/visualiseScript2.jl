# using SparseArrays
using UnPack
using CairoMakie 
using FromFile
using DrWatson
using OrdinaryDiffEq
# using Printf
# using Dates
# using InvertedIndices
# using Statistics
# using JLD2

@from "$(srcdir("Visualise.jl"))" using Visualise

subFolder = "CSF"
folderName = "24-12-18-15-07-24_K₂=0.12_K₄=0.1_Tᵣ=0.005_differencing=centre_nSpatialDims=1_thicknessProfile=GRF_α_C=5.0_β=8.8_𝓓=34.6"
# Create frames subdirectory to store system state at each output time

@unpack rawParams, derivedParams = (load(datadir("sims", subFolder, folderName, "params.jld2")))

files = [f for f in readdir(datadir("sims", subFolder, folderName)) if occursin("results", f)]
u = []
for file in files
    data = load(datadir("sims", subFolder, folderName, file))
    push!(u, data["u"])
end

# @unpack p = load(datadir("sims", subFolder, folderName, files[1]))

#%%

if length(p.dims)==2
    concentrationSurfaceMovie(u, p.dims; subFolder=subFolder, folderName=folderName)
    concentrationHeatmapMovie(u, p.dims; subFolder=subFolder, folderName=folderName)
    spaceIntegralOver_ν_Movie(u, p; subFolder=subFolder, folderName=folderName)
    # if thicknessProfile=="GRF"
        # thicknessPlot(p.hᵥ, p.dims; subFolder=subFolder, folderName=folderName)
    # end
else
    spaceIntegralOver_ν_Movie(u, p; subFolder=subFolder, folderName=folderName)
    uSlices = [selectdim(reshape(u, p.dims...), 3, dims[3]÷2) for u in u]
    uSlicesReshaped = [reshape(u, prod(dims[Not(3)])) for u in uSlices]
    concentrationSurfaceMovie(uSlicesReshaped, p.dims; subFolder=subFolder, folderName=folderName)
    concentrationHeatmapMovie(uSlicesReshaped, p.dims; subFolder=subFolder, folderName=folderName)
    # if thicknessProfile=="GRF"
        # thicknessPlot(p.hᵥ, p.dims; subFolder=subFolder, folderName=folderName)
    # end
end