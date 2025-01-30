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

textheme = Theme(fonts=(; regular=texfont(:text),
                        bold=texfont(:bold),
                        italic=texfont(:italic),
                        bold_italic=texfont(:bolditalic)))

@from "$(srcdir("Glycosylation.jl"))" using Glycosylation
@from "$(srcdir("Visualise.jl"))" using Visualise
@from "$(srcdir("UsefulFunctions.jl"))" using UsefulFunctions
@from "$(srcdir("MakeWeightMatrices.jl"))" using MakeWeightMatrices
@from "$(srcdir("DerivedParameters.jl"))" using DerivedParameters

#%%

subFolder = "analyticNumericFit2"
terminateAt = "nuWall"
thicknessProfile = "uniform"
differencing = "centre"
solver = SSPRK432()
nOutputs = 1000
σGRF = 0.2

nSpatialDims = 1
Ngrid = 401
dims = fill(Ngrid, nSpatialDims+1)

include(projectdir("notebooks", "paramsRaw.jl"))

#%%

derivedParams = derivedParameters(Ω, Ωperp, N, k_Cd, k_Ca, k_Sd, k_Sa, k₁, k₂, k₃, k₄, 𝒞, 𝒮, ℰ, D_C, D_S, Tᵣstar; checks=true)
@unpack L₀, E₀, C_b, S_b, δ_C, δ_S, α_C, α_S, C₀, S₀, Tᵣ, T̃ᵣ, K₂, K₃, K₄, σ, ϵ, 𝒟, β, h_C, h_S,u, λ, ζ, γ, Δ, F = derivedParams

#%%

# Create directory for run data labelled with current time.
paramsName = @savename nSpatialDims K₂ K₄ α_C β 𝒟 T̃ᵣ thicknessProfile differencing
folderName = "$(Dates.format(Dates.now(),"yy-mm-dd-HH-MM-SS"))_$(paramsName)"
# Create frames subdirectory to store system state at each output time
mkpath(datadir("sims",subFolder,folderName))

#%%

sol, p = glycosylationAnyD(dims, K₂, K₄, T̃ᵣ, α_C, 𝒟, β, thickness=thicknessProfile, differencing=differencing, solver=solver, nOutputs=nOutputs, terminateAt=terminateAt)
println("finished sim")

#%%

νMax = 1.0
νs   = collect(range(0.0, νMax, dims[1]))

#%%

midpoint = length(sol.u)÷2
C_peak, ind_peak = findmax(reshape(sol.u[midpoint], dims...)[:,1])
ν_peak = νs[ind_peak]
Ẽ = K₂/(1+K₂)
D = Ẽ*K₂*K₄/(1+α_C)
t₀ = sol.t[midpoint] - 1/(4.0*π*D*C_peak^2)
ν₀ = ν_peak - Ẽ*β*(sol.t[midpoint]-t₀)/(1+α_C)
νsOffset = νs.-ν₀
tsOffset = sol.t.-t₀
firstPositivetIndex = findfirst(x->x>0, tsOffset)

#%%

fig = Figure(size=(1000,1000), fontsize=32, theme=textheme)

ax1 = CairoMakie.Axis(fig[1, 1])
ax1.xlabel = L"\nu"
ax1.ylabel = L"\tilde{C}"
analyticLine = Observable(zeros(dims[1]))
numericLine = Observable(zeros(dims[1]))
l1 = lines!(ax1, νs, analyticLine, color=:red, linewidth=4)
l2 = lines!(ax1, νs, numericLine, color=:blue, linewidth=4)
Legend(fig[1,2], [l1, l2], ["Analytic", "Numeric"])
ylims!(ax1, (-2.0, maximum(sol.u[1])))
xlims!(ax1, (0.0, 1.0))
analyticVals = zeros(size(νsOffset)) 

ax2 = CairoMakie.Axis(fig[2, 1])#, aspect=1)
ax2.xlabel = L"\tilde{t}"
ax2.ylabel = L"\tilde{M}_\phi"
xlims!(ax2, (0.0, sol.t[end]))
ylims!(ax2, (0.0, 1.05*π))
Ms = Observable(zeros(length(sol.t)))
Ts = Observable(zeros(length(sol.t)))
l3 = lines!(ax2, Ts, Ms, linewidth=4)

record(fig, datadir("sims",subFolder, folderName, "analyticCs.mp4"), 1:length(sol.t); framerate=20) do i
    if tsOffset[i] > 0
        analyticVals .= homogeneousWidthC.(νsOffset, K₂, K₄, α_C, β, tsOffset[i])
        analyticLine[] .= analyticVals
        uInternal = reshape(sol.u[i], dims...)
        numericLine[] .= uInternal[:,dims[2]÷2]
        analyticLine[] = analyticLine[]
        numericLine[] = numericLine[]

        Ts[][i] = sol.t[i]
        Ts[] = Ts[]   
        Mϕ = M̃ϕ(sol.u[i], p.W, p.dims, p.dν, p.hᵥ, ϕ)
        Ms[][i] = Mϕ
        Ms[] = Ms[]
    end
end

#%%

fig = Figure(size=(1000,1000), fontsize=32, theme=textheme)
ax1 = CairoMakie.Axis(fig[1, 1])

allLines = []
allTs = []
colorsUsed = [(:red), (:green), (:blue)]
for (c,i) in enumerate([firstPositivetIndex, (length(sol.t)-firstPositivetIndex)÷2+firstPositivetIndex, length(sol.t)])
    uInternal = reshape(sol.u[i], dims...)
    push!(allLines, lines!(ax1, νs, uInternal[:,1], linestyle=:solid, color=(colorsUsed[c], 0.5), linewidth=4))
    push!(allLines, lines!(ax1, νs, homogeneousWidthC.(νsOffset, K₂, K₄, α_C, β, tsOffset[i]), linestyle=:dot, color=(colorsUsed[c], 1.0), linewidth=4))
    push!(allTs, @sprintf("%.1f", tsOffset[i]))
end

labels = []
for t in allTs
    push!(labels, "Numeric, t=$t")
    push!(labels, "Analytic, t=$t")
end

axislegend(ax1, allLines, labels, labelsize = 16)
ax1.xlabel = L"\nu"
ax1.ylabel = L"\tilde{C}"
ylims!(ax1, (0.0, 20.0))#maximum(sol.u[1])))
xlims!(ax1, (0.0, 1.0))

ax2 = Axis(fig[2,1], yticks = (0.0:π/2.0:π, [L"0", L"π/2", L"π"]))
ylims!(ax2, (0.0, 1.05*π))
xlims!(ax2, (0.0, sol.t[end]))
l3 = lines!(ax2, Ts, Ms[], linewidth=4, color=(:black, 1.0))

ax2.xlabel = L"\tilde{t}"
ax2.ylabel = L"\tilde{M}_\phi"

save(datadir("sims", subFolder, folderName, "analyticComparisonν_0=$(@sprintf("%f", ν₀))t_0=$(@sprintf("%f", t₀)).png"), fig)
display(fig)

@show t₀
@show ν₀

