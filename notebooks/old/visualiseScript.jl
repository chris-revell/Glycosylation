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
folderName = "24-12-16-18-31-33_K₂=0.12_K₄=0.1_Tᵣ=5.0_differencing=centre_nSpatialDims=1_thicknessProfile=GRF_α_C=5.0_β=8.8_𝓓=34.6"
# Create frames subdirectory to store system state at each output time
data = load(datadir("sims", subFolder, folderName, "output.jld2"))
@unpack sol, p = data

#%%

if length(p.dims)==2
    concentrationSurfaceMovie(sol.u, p.dims; subFolder=subFolder, folderName=folderName)
    concentrationHeatmapMovie(sol.u, p.dims; subFolder=subFolder, folderName=folderName)
    M̃movie(sol.u, p; subFolder=subFolder, folderName=folderName)
    # if thicknessProfile=="GRF"
        # thicknessPlot(p.hᵥ, p.dims; subFolder=subFolder, folderName=folderName)
    # end
else
    M̃movie(sol.u, p; subFolder=subFolder, folderName=folderName)
    uSlices = [selectdim(reshape(u, p.dims...), 3, dims[3]÷2) for u in sol.u]
    uSlicesReshaped = [reshape(u, prod(dims[Not(3)])) for u in uSlices]
    concentrationSurfaceMovie(uSlicesReshaped, p.dims; subFolder=subFolder, folderName=folderName)
    concentrationHeatmapMovie(uSlicesReshaped, p.dims; subFolder=subFolder, folderName=folderName)
    # if thicknessProfile=="GRF"
        # thicknessPlot(p.hᵥ, p.dims; subFolder=subFolder, folderName=folderName)
    # end
end