
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
@from "$(srcdir("CisternaWidth.jl"))" using CisternaWidth

#%%

subFolder = ""

thicknessProfile = "uniform"
differencing = "centre"
solver = SSPRK432()
nOutputs = 100
σGRF = 0.2

nSpatialDims = 1
Ngrid = 401
dims = fill(Ngrid, nSpatialDims+1)

#%%

h₀ = 0.002
Ωperp = 10000    # Dimensional lumen footprint area
Ω     = h₀*Ωperp      # Dimensional lumen volume 
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
Tᵣstar= 10000000000000.0  # Release time
ϕ     = 0.5

#%%

derivedParams = derivedParameters(Ω, Ωperp, N, k_Cd, k_Ca, k_Sd, k_Sa, k₁, k₂, k₃, k₄, 𝓒, 𝓢, 𝓔, D_C, D_S, Tᵣstar; checks=true)
@unpack L₀, E₀, C_b, S_b, δ_C, δ_S, α_C, α_S, C₀, S₀, Tᵣ, T̃ᵣ, K₂, K₃, K₄, σ, ϵ, 𝓓, β = derivedParams

#%%

# Create directory for run data labelled with current time.
paramsName = @savename nSpatialDims K₂ K₄ α_C β 𝓓 Tᵣ thicknessProfile differencing
folderName = "$(Dates.format(Dates.now(),"yy-mm-dd-HH-MM-SS"))_$(paramsName)"
# Create frames subdirectory to store system state at each output time
mkpath(datadir("sims",subFolder,folderName))

#%%

sol, p = glycosylationAnyD(dims, K₂, K₄, T̃ᵣ, α_C, 𝓓, β, thickness=thicknessProfile, differencing=differencing, solver=solver, nOutputs=nOutputs, σGRF=σGRF)
println("finished sim")

#%%

rawParams = (
    thicknessProfile = thicknessProfile,
    differencing = differencing,
    solver = solver,
    nOutputs = nOutputs,
    σGRF = σGRF,
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
    𝓒 = 𝓒,
    𝓢 = 𝓢,
    𝓔 = 𝓔,
    D_C = D_C,
    D_S = D_S,
    Tᵣstar = Tᵣstar,
    ϕ = ϕ
)

jldsave(datadir("sims",subFolder,folderName,"solution.jld2"); sol, p, rawParams)

#%%

if nSpatialDims==1
    concentrationSurfaceMovie(sol.u, dims; subFolder=subFolder, folderName=folderName)
    concentrationHeatmapMovie(sol.u, dims; subFolder=subFolder, folderName=folderName)
    spaceIntegralOver_ν_Movie(sol.u, p; subFolder=subFolder, folderName=folderName)
    if thicknessProfile=="GRF"
        thicknessPlot(p.hᵥ, p.dims; subFolder=subFolder, folderName=folderName)
    end
else
    spaceIntegralOver_ν_Movie(sol.u, p; subFolder=subFolder, folderName=folderName)
    uSlices = [selectdim(reshape(u, dims...), 3, dims[3]÷2) for u in sol.u]
    uSlicesReshaped = [reshape(u, prod(dims[Not(3)])) for u in uSlices]
    concentrationSurfaceMovie(uSlicesReshaped, dims; subFolder=subFolder, folderName=folderName)
    concentrationHeatmapMovie(uSlicesReshaped, dims; subFolder=subFolder, folderName=folderName)
    if thicknessProfile=="GRF"
        thicknessPlot(hᵥ, dims; subFolder=subFolder, folderName=folderName)
    end
end


# fig = Figure(size=(500,500))

# ax2 = CairoMakie.Axis(fig[1, 1])#, aspect=1)
# ax2.xlabel = "t"
# ax2.ylabel = L"M_\phi"
# xlims!(ax2, (0.0, sol.t[end]))
# ylims!(ax2, (0.0, M_star_ϕ(sol.u[end], p.W, dims, p.dν, p.hᵥ, α_C, C_b, Ω, ϕ)))

# Ms = Observable(zeros(length(sol.t)))
# Ts = Observable(zeros(length(sol.t)))
# l3 = lines!(ax2, Ts, Ms, linewidth=4)

# record(fig, datadir("sims",subFolder,folderName,"usefulProduction.mp4"), 1:length(sol.t); framerate=50) do i
#     Ts[][i] = sol.t[i]
#     Ts[] = Ts[]        
#     Mϕ = M_star_ϕ(sol.u[i], p.W, dims, p.dν, p.hᵥ, α_C, C_b, Ω, ϕ)
#     Ms[][i] = Mϕ
#     Ms[] = Ms[]
# end


# W = vertexVolumeWeightsMatrix(dims, spacing)
# hᵥ = spdiagm(reshape(mat_h, prod(dims)))
# h₀ = 0.01
# Ωperp = 10000
# Ω = Ωperp*h₀
# N = 100
# 𝓒     = 100000.0
# 𝓢     = 100000.0
# 𝓔     = 0.0001
# C_b  = 𝓒/Ω
# k₁ = 1.0
# @show P_star(sol.u[end], W, dims, dν, hᵥ, α_C, C_b, Ω, ϕ, Ωperp, k₁, 𝓔, Tᵣ)
# @show 𝓟starUniform(𝓒, 𝓔, 𝓢, ϕ, N, k₁, k₂, K₃, K₄, Ωperp, h₀, h_C, h_S)

# Mϕ50 = α_C*𝓒/(2.0*(1+α_C))


# T50 = 2.0*Ωperp/(k₁*𝓔) * N^2* (K₂+σ*K₃) * (t_0 + (ϕ-ν_0)*(1+α_C)*())
# P50 = Mϕ50/T50


# Tᵣ = 30.0
# K₂ = 1.0
# K₄ = 0.0001
# α_C = 1.0
# 𝓓 = 1.0
# β = 0.1

# differencing = "upstream"

# include(projectdir("notebooks","paramsRaw.jl"))
# h_C = 2*k_Ca/k_Cd
# h_S = 2*k_Sa/k_Sd
# hMax = h_C*25
# hMin = h_C/10
# h₀s = collect(hMin:2*hMin:hMax)
# Ωs    = h₀s.*Ωperp      # Dimensional lumen volume 

# h₀ = 2*h_C
# Ω = h₀*Ωperp
# nSpatialDims == 1 ? spacing  = [dν, dx] : spacing  = [dν, dx, dy]