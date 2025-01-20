
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
nOutputs = 500
σGRF = 0.2

nSpatialDims = 1
Ngrid = 401
dims = fill(Ngrid, nSpatialDims+1)

include(projectdir("notebooks", "paramsRaw.jl"))


#%%

derivedParams = derivedParameters(Ω, Ωperp, N, k_Cd, k_Ca, k_Sd, k_Sa, k₁, k₂, k₃, k₄, 𝒞, 𝒮, ℰ, D_C, D_S, Tᵣstar; checks=true)
@unpack L₀, E₀, C_b, S_b, δ_C, δ_S, α_C, α_S, C₀, S₀, Tᵣ, T̃ᵣ, K₂, K₃, K₄, σ, ϵ, 𝒟, β, h_C, h_S = derivedParams

#%%

sol, p = glycosylationAnyD(dims, K₂, K₄, T̃ᵣ, α_C, 𝒟, β, thickness=thicknessProfile, differencing=differencing, solver=solver, nOutputs=nOutputs, terminateAt=terminateAt)
println("finished sim")

#%%

# Create directory for run data labelled with current time.
paramsName = @savename nSpatialDims K₂ K₄ α_C β 𝒟 T̃ᵣ thicknessProfile differencing
folderName = "$(Dates.format(Dates.now(),"yy-mm-dd-HH-MM-SS"))_$(paramsName)"
# Create frames subdirectory to store system state at each output time
mkpath(datadir("sims",subFolder,folderName))

#%%

# xMax = π^(1/nSpatialDims)
# xs   = collect(range(0.0, xMax, dims[2]))
# dx   = xs[2]-xs[1]
# if nSpatialDims > 1 
#     yMax = xMax
#     ys   = collect(range(0.0, yMax, dims[3]))
#     dy   = ys[2]-ys[1]
# end
νMax = 1.0
νs   = collect(range(0.0, νMax, dims[1]))
# dν   = νs[2]-νs[1]
# nSpatialDims == 1 ? spacing  = [p.dν, dx] : spacing  = [p.dν, dx, dy]

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

fig = Figure(size=(1000,1000))
ax = CairoMakie.Axis(fig[1, 1])#, aspect=1)
ax.xlabel = "ν"
ax.ylabel = "C"
analyticLine = Observable(zeros(dims[1]))
numericLine = Observable(zeros(dims[1]))
l1 = lines!(ax, νs, analyticLine, color=:red)
l2 = lines!(ax, νs, numericLine, color=:blue)
Legend(fig[1,2], [l1, l2], ["Analytic", "Numeric"])
ylims!(ax, (-2.0, 20.0))
xlims!(ax, (0.0, 1.0))
analyticVals = zeros(size(νsOffset)) # homogeneousWidthC.(νsOffset, K₂, K₄, α_C, β, tsOffset[1])

ax2 = CairoMakie.Axis(fig[2, 1])#, aspect=1)
ax2.xlabel = "t"
ax2.ylabel = "Mϕ"
xlims!(ax2, (0.0, sol.t[end]))
ylims!(ax2, (0.0, M_star_ϕ(sol.u[end], p.W, dims, p.dν, p.hᵥ, α_C, C_b, Ω, ϕ)))
Ms = Observable(zeros(length(sol.t)))
Ts = Observable(zeros(length(sol.t)))
l3 = lines!(ax2, Ts, Ms)

record(fig, datadir("sims",subFolder, folderName, "analyticCs.mp4"), 1:length(sol.t); framerate=50) do i
    if tsOffset[i] > 0
        analyticVals .= homogeneousWidthC.(νsOffset, K₂, K₄, α_C, β, tsOffset[i])
        analyticLine[] .= analyticVals
        uInternal = reshape(sol.u[i], dims...)
        numericLine[] .= uInternal[:,dims[2]÷2]
        analyticLine[] = analyticLine[]
        numericLine[] = numericLine[]

        Ts[][i] = sol.t[i]
        Ts[] = Ts[]        
        Mϕ = M_star_ϕ(sol.u[i], p.W, dims, p.dν, p.hᵥ, α_C, C_b, Ω, ϕ)
        Ms[][i] = Mϕ
        Ms[] = Ms[]
    end
end

#%%



#%%

fig = Figure(size=(1000,1000), fontsize=32)
ax = CairoMakie.Axis(fig[1, 1])
# ylims!(ax, (-2.0, 20.0))
ylims!(ax, (0.0, 20.0))
xlims!(ax, (0.0, 1.0))

allLines = []
analyticLines = []
allTs = []

colorsUsed = [(:red), (:green), (:blue)]

for (c,i) in enumerate([firstPositivetIndex, (length(sol.t)-firstPositivetIndex)÷2+firstPositivetIndex, length(sol.t)])
    uInternal = reshape(sol.u[i], dims...)
    push!(allLines, lines!(ax, νs, uInternal[:,1], linestyle=:solid, color=(colorsUsed[c], 0.5), linewidth=4))
    push!(allLines, lines!(ax, νs, homogeneousWidthC.(νsOffset, K₂, K₄, α_C, β, tsOffset[i]), linestyle=:dot, color=(colorsUsed[c], 1.0), linewidth=4))
    push!(allTs, @sprintf("%.1f", tsOffset[i]))
end

labels = []
for t in allTs
    push!(labels, "Numeric, t=$t")
    push!(labels, "Analytic, t=$t")
end

# Legend(fig[1,1], allLines, labels, labelsize = 2)
axislegend(ax, allLines, labels, labelsize = 16)
ax.xlabel = L"\nu"
ax.ylabel = L"C"

ax2 = Axis(fig[2,1])
# ylims!(ax2, (0.0, M_star_ϕ(sol.u[end], p.W, dims, p.dν, p.hᵥ, α_C, C_b, Ω, ϕ)))
l3 = lines!(ax2, Ts, Ms[]./Ms[][end], linewidth=4, color=(:black, 1.0))
ax2.xlabel = "Time"
ax2.ylabel = L"M^*_\phi"

save(datadir("sims", subFolder, folderName, "analyticComparisonν_0=$(@sprintf("%f", ν₀))t_0=$(@sprintf("%f", t₀)).png"), fig)
display(fig)

@show β
@show α_C
@show t₀
@show ν₀
@show T̃ᵣ
@show K₂
@show K₄
@show 𝒟
@show ϕ


# Tᵣ = 30.0
# K₂ = 1.0
# K₄ = 0.0001
# α_C = 1.0
# 𝒟 = 1.0
# β = 0.1

# #%%
# thicknessProfile = "uniform"
# differencing = "centre"
# nSpatialDims = 1
# Ngrid = 401
# # dims = [Ngrid,2]
# dims = fill(Ngrid, nSpatialDims+1)

# #%%

# h₀ = 0.1
# Ωperp = 10000    # Dimensional lumen footprint area
# Ω     = h₀*Ωperp      # Dimensional lumen volume 
# N     = 100     # Maximum polymer length 
# k_Cd  = 1.0 # Complex desorption rate
# k_Ca  = 0.01 # Complex adsorption rate
# k_Sd  = 1.0 # Substrate desorption rate
# k_Sa  = 0.01 # Substrate adsorption rate
# k₁    = 1.0   # Complex formation forward reaction rate 
# k₂    = 0.1   # Complex dissociation reverse reaction rate 
# k₃    = 0.1   # Product formation
# k₄    = 0.1  # Product dissociation 
# 𝒞     = 100000.0
# 𝒮     = 100000.0
# ℰ     = 0.0001
# D_C   = 0.0000001  # Monomer/polymer diffusivity
# D_S   = 0.0000001  # Substrate diffusivity
# Tᵣstar= 1000000000.0  # Release time
# ϕ     = 0.5
