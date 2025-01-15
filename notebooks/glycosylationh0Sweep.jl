
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

thicknessProfile = "Gaussian"
differencing = "centre"
solver = SSPRK432()
nOutputs = 100

Ωperp = 10000    # Dimensional lumen footprint area
N     = 100     # Maximum polymer length 
k_Cd  = 1.0 # Complex desorption rate
k_Ca  = 0.01 # Complex adsorption rate
k_Sd  = 1.0 # Substrate desorption rate
k_Sa  = 0.01 # Substrate adsorption rate
k₁    = 1.0   # Complex formation forward reaction rate 
k₂    = 0.1   # Complex dissociation reverse reaction rate 
k₃    = 0.1   # Product formation
k₄    = 0.1  # Product dissociation 
𝓒     = 100000.0
𝓢     = 100000.0
𝓔     = 0.0001
D_C   = 0.0000001  # Monomer/polymer diffusivity
D_S   = 0.0000001  # Substrate diffusivity
Tᵣstar= 100000000000.0  # Release time
ϕ     = 0.5

nSpatialDims = 1
Ngrid = 401
# dims = fill(Ngrid, nSpatialDims+1)
dims = [Ngrid,2]

include(projectdir("notebooks","paramsRaw.jl"))
h_C = 2*k_Ca/k_Cd
h_S = 2*k_Sa/k_Sd
hMax = h_C*15
hMin = h_C/10
h₀s = collect(hMin:hMin:hMax)
Ωs    = h₀s.*Ωperp      # Dimensional lumen volume 

#%%

# include(projectdir("notebooks","paramsDerived.jl"))

#%%

xMax = π^(1/nSpatialDims)
xs   = collect(range(0.0, xMax, dims[2]))
dx   = xs[2]-xs[1]
if nSpatialDims > 1 
    yMax = xMax
    ys   = collect(range(0.0, yMax, dims[3]))
    dy   = ys[2]-ys[1]
end
νMax = 1.0
νs   = collect(range(0.0, νMax, dims[1])) # Positions of discretised vertices in polymerisation space 
dν   = νs[2]-νs[1]
nSpatialDims == 1 ? spacing  = [dν, dx] : spacing  = [dν, dx, dy]

# W = vertexVolumeWeightsMatrix(dims, spacing)

PstarsAnalytic = []
PstarsSim = []
MstarsPhiSim = []
sols = []
for i=1:length(h₀s)
    @show h₀s[i]
    
    derivedParams = derivedParameters(Ωs[i], Ωperp, N, k_Cd, k_Ca, k_Sd, k_Sa, k₁, k₂, k₃, k₄, 𝓒, 𝓢, 𝓔, D_C, D_S, Tᵣstar; checks=false)
    @unpack L₀, E₀, C_b, S_b, δ_C, δ_S, α_C, α_S, C₀, S₀, Tᵣ, T̃ᵣ, K₂, K₃, K₄, σ, ϵ, 𝓓, β = derivedParams

    sol, p = glycosylationAnyD(dims, K₂, K₄, T̃ᵣ, α_C, 𝓓, β, thickness=thicknessProfile, differencing=differencing, solver=solver, nOutputs=nOutputs, terminateAt="halfProduction")

    # push!(PstarsAnalytic, 𝓟starUniform(𝓒, 𝓔, 𝓢, ϕ, N, k₁, k₂, K₃, K₄, Ωperp, h₀s[i], h_C, h_S))
end

#%%

linesVec = []
labelsVec = []
fig = Figure()#size=(500,500))
# ax1 = Axis(fig[1,1])
# push!(linesVec, lines!(ax1, h₀s, MstarsPhiSim, color=:blue))
# push!(labelsVec, "Numerical")
# ylims!(ax1, (0.0,maximum(MstarsPhiSim)))
# xlims!(ax1, (0.0,maximum(h₀s)))
# ax1.xlabel = "h₀"
# ax1.ylabel = L"M^*"

ax2 = Axis(fig[1,1])
push!(linesVec, lines!(ax2, h₀s, PstarsSim, color=:blue))
push!(labelsVec, "Numerical")
ylims!(ax2, (0.0,maximum(PstarsSim)))
xlims!(ax2, (0.0,maximum(h₀s)))
ax2.xlabel = "h₀"
ax2.ylabel = L"𝓟^*"

# ax3 = Axis(fig[3,1])
# ax3.xlabel = "h₀"
# ax3.ylabel = L"𝓟^*"
# ylims!(ax3, (0.0,maximum(PstarsAnalytic)))
# xlims!(ax3, (0.0,maximum(h₀s)))
# push!(linesVec, lines!(ax3, h₀s, PstarsAnalytic, color=:red))
# push!(labelsVec, "Analytic")

# push!(linesVec, vlines!(ax1, h_C, color=:green))
# push!(labelsVec, L"h_C")
# push!(linesVec, vlines!(ax2, h_C, color=:green))
# push!(labelsVec, L"h_C")
# push!(linesVec, vlines!(ax1, h_S, color=:orange))
# push!(labelsVec, L"h_S")
# push!(linesVec, vlines!(ax2, h_S, color=:orange))
# push!(labelsVec, L"h_S")

# Legend(fig[:,2], linesVec, labelsVec)

display(fig)

paramsName = @savename k_Cd k_Ca k_Sd k_Sa k₁ k₂ k₃ k₄ 𝓒 𝓢 𝓔 D_C D_S Tᵣstar ϕ
folderName = "$(Dates.format(Dates.now(),"yy-mm-dd-HH-MM-SS"))_$(paramsName)"
save("$(Dates.format(Dates.now(),"yy-mm-dd-HH-MM-SS"))_$(paramsName)_simulationPvsh.png",fig)



