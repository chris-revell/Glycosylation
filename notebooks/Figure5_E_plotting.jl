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


νs = collect(range(0.0,1.0,p1.dims[1]))
xs = collect(range(0.0,sqrt(π),p1.dims[2]))
ys = collect(range(0.0,sqrt(π),p1.dims[2]))

#%%

Ẽs = []
for u in sol1.u 
    Ẽ = zeros(size(p1.matE))
    uMat = reshape(u, p1.dims...)
    integ = p1.dν.*(sum(uMat, dims=1) .- 0.5.*selectdim(uMat, 1, 1) .- 0.5.*selectdim(uMat, 1, 1))
    for slice in eachslice(Ẽ, dims=1)
        slice .= p1.matFₑ.*(p1.K₂./(p1.K₂ .+ selectdim(integ, 1, 1)))
    end
    push!(Ẽs, Ẽ[1,:,:])
end

minimum(Ẽs[27])

fig = Figure(size=(500,500))
ax = Axis(fig[1, 1], aspect=DataAspect())
ax.xlabel = L"x"
ax.ylabel = L"y"
plottedẼ = Observable(copy(Ẽs[27]))
# globalmin = minimum([minimum(Ẽ) for Ẽ in Ẽs])
globalmin = minimum(Ẽs[27])
# globalmax = maximum([maximum(Ẽ) for Ẽ in Ẽs])
globalmax = maximum(Ẽs[27])
clims = (globalmin,globalmax)
xs = collect(range(0,sqrt(π),dims[2]))
ys = collect(range(0,sqrt(π),dims[2]))
heatmap!(ax, xs, ys, plottedẼ, colorrange=clims, colormap=:batlow)
Colorbar(fig[1,2], limits=clims, label=L"\tilde{E}")

display(fig)


record(fig, datadir("sims",subFolder,folderName,"E_HeatmapMovie.mp4"), 1:length(Ẽs); framerate=10) do i
    plottedẼ[] .= Ẽs[i]
    plottedẼ[] = plottedẼ[]
end


#%%

testFrame = findfirst(x->x>0.1, sol1.t)

maximum(Ẽs[testFrame])
mean(Ẽs[testFrame])
minimum(Ẽs[testFrame])

100.0*minimum(Ẽs[testFrame])/mean(Ẽs[testFrame])
100.0*maximum(Ẽs[testFrame])/mean(Ẽs[testFrame])