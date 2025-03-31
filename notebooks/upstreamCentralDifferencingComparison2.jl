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

subFolder = "upstreamVsCentralDifferencing"
terminateAt = "nuWall"
thicknessProfile = "uniform"
# differencing = "centre"
solver = SSPRK432()
nOutputs = 200

nSpatialDims = 1
Ngrid = 201
dims = fill(Ngrid, nSpatialDims+1)

include(projectdir("notebooks", "paramsRaw.jl"))


derivedParams = derivedParameters(Ω, Ωperp, N, k_Cd, k_Ca, k_Sd, k_Sa, k₁, k₂, k₃, k₄, 𝒞, 𝒮, ℰ, D_C, D_S, Tᵣstar; checks=true)
@unpack L₀, E₀, C_b, S_b, δ_C, δ_S, α_C, α_S, C₀, S₀, Tᵣ, T̃ᵣ, K₂, K₃, K₄, σ, ϵ, 𝒟, β, h_C, h_S, λ, ζ, γ, Δ, F = derivedParams

#%%

# Create directory for run data labelled with current time.
paramsName = @savename nSpatialDims K₂ K₄ α_C β 𝒟 T̃ᵣ thicknessProfile 
folderName = "$(Dates.format(Dates.now(),"yy-mm-dd-HH-MM-SS"))_$(paramsName)"
# Create frames subdirectory to store system state at each output time
mkpath(datadir("sims",subFolder,folderName))

#%%

solCentre, p = glycosylationAnyD(dims, K₂, K₄, T̃ᵣ, α_C, 𝒟, β, thickness="uniform", differencing="centre", solver=solver, nOutputs=nOutputs, terminateAt=terminateAt)
solUpstream, p = glycosylationAnyD(dims, K₂, K₄, T̃ᵣ, α_C, 𝒟, β, thickness="uniform", differencing="upstream", solver=solver, nOutputs=nOutputs, terminateAt=terminateAt)
println("finished sim")

#%%

rawParams = (
    thicknessProfile = thicknessProfile,
    # differencing = differencing,
    solver = solver,
    nOutputs = nOutputs,
    # σGRF = σGRF,
    # λGRF = λGRF,
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

jldsave(datadir("sims",subFolder,folderName,"solution.jld2"); solCentre, solUpstream, p, rawParams)


νs   = collect(range(0.0, 1.0, dims[1])) 

fig = Figure(size=(1000,1000))
ax = CairoMakie.Axis(fig[1, 1], aspect=1)
ax.xlabel = "ν"
ax.ylabel = "C"
analyticLine = Observable(zeros(dims[1]))
numericLine = Observable(zeros(dims[1]))
l1 = lines!(ax, νs, analyticLine, color=(:red, 0.5), linewidth=4, linestyle=:dot)
l2 = lines!(ax, νs, numericLine, color=(:blue, 0.5), linewidth=4)#, linestyle=:dot)
Legend(fig[1,2], [l1, l2], ["Upstream differencing", "Central differencing"])
ylims!(ax, (-20.0, 20.0))
xlims!(ax, (0.0, 1.0))
# analyticVals = homogeneousWidthC.(νsOffset, K₂, K₄, α_C, β, tsOffset[1])
record(fig, datadir("sims",subFolder, folderName, "analyticCs.mp4"), 1:length(solUpstream.t); framerate=50) do i
    uInternal = reshape(solUpstream.u[i], dims...)
    analyticLine[] .= uInternal[:,1]
    uInternal = reshape(solCentre.u[i], dims...)
    numericLine[] .= uInternal[:,1]
    analyticLine[] = analyticLine[]
    numericLine[] = numericLine[]
    # if i in [1, 251, 500]
    #     save(datadir("sims",subFolder, folderName, "analyticCs$i.png"), fig)
    # end
end

#%%

midpoint = length(solCentre.u)÷2
C_peak, ind_peak = findmax(reshape(solCentre.u[midpoint], dims...)[:,1])
ν_peak = νs[ind_peak]
Ẽ = K₂/(1+K₂)
D = Ẽ*K₂*K₄/(1+α_C)
t₀Centre = solCentre.t[midpoint] - 1/(4.0*π*D*C_peak^2)
ν₀Centre = ν_peak - Ẽ*β*(solCentre.t[midpoint]-t₀Centre)/(1+α_C)

νsOffsetCentre = νs.-ν₀Centre
tsOffsetCentre = solCentre.t.-t₀Centre

midpoint = length(solUpstream.u)÷2
C_peak, ind_peak = findmax(reshape(solUpstream.u[midpoint], dims...)[:,1])
ν_peak = νs[ind_peak]
Ẽ = K₂/(1+K₂)
D = Ẽ*K₂*K₄/(1+α_C)
t₀Upstream = solUpstream.t[midpoint] - 1/(4.0*π*D*C_peak^2)
ν₀Upstream = ν_peak - Ẽ*β*(solUpstream.t[midpoint]-t₀Upstream)/(1+α_C)

νsOffsetUpstream = νs.-ν₀Upstream
tsOffsetUpstream = solUpstream.t.-t₀Upstream

#%%

fig = Figure(size=(1000,750), fontsize=32)
ax = CairoMakie.Axis(fig[1, 1])
ylims!(ax, (0.0, 20.0))
xlims!(ax, (0.0, 1.0))
allLines = []
for (c,i) in enumerate([1, 251, 500])
    uInternal1 = reshape(solUpstream.u[i], dims...)
    push!(allLines, lines!(ax, νs, uInternal1[:,dims[2]÷2], linestyle=:solid, color=(:blue, 1.0), linewidth=4))
    uInternal2 = reshape(solCentral.u[i], dims...)
    push!(allLines, lines!(ax, νs, uInternal2[:,dims[2]÷2], linestyle=:solid, color=(:blue, 1.0), linewidth=4))
    push!(allLines, lines!(ax, νs, homogeneousWidthC.(νsOffset, K₂, K₄, α_C, β, tsOffset[i]), linestyle=:dash, color=(:red,1.0), linewidth=4))
end
Legend(fig[1,2], allLines[1:2], ["Numeric", "Analytic"])
ax.xlabel = L"\nu"
ax.ylabel = L"C"
save(datadir("sims", subFolder, folderName, "analyticComparison.png"), fig)

# #%%

# fig = Figure(size=(1000,1000))
# ax = CairoMakie.Axis(fig[1, 1], aspect=1)
# ax.xlabel = "ν"
# ax.ylabel = "C"
# analyticLine = Observable(zeros(dims[1]))
# numericLine = Observable(zeros(dims[1]))
# l1 = lines!(ax, νs, analyticLine, color=:red)
# l2 = lines!(ax, νs, numericLine, color=:blue)
# Legend(fig[1,2], [l1, l2], ["Analytic", "Numeric"])
# ylims!(ax, (-20.0, 20.0))
# xlims!(ax, (0.0, 1.0))
# analyticVals = homogeneousWidthC.(νsOffset, K₂, K₄, α_C, β, tsOffset[1])
# record(fig, datadir("sims",subFolder, folderName, "analyticCs.mp4"), 1:length(sol.t); framerate=50) do i
#     analyticVals .= homogeneousWidthC.(νsOffset, K₂, K₄, α_C, β, tsOffset[i])
#     analyticLine[] .= analyticVals
#     uInternal = reshape(sol.u[i], dims...)
#     numericLine[] .= uInternal[:,dims[2]÷2]
#     analyticLine[] = analyticLine[]
#     numericLine[] = numericLine[]
#     # if i in [1, 251, 500]
#     #     save(datadir("sims",subFolder, folderName, "analyticCs$i.png"), fig)
#     # end
# end

#%%