
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
terminateAt = "halfProduction"
thicknessProfile = "uniform"
differencing = "centre"
solver = SSPRK432()
nOutputs = 100
σGRF = 0.2

nSpatialDims = 2
Ngrid = 101
dims = fill(Ngrid, nSpatialDims+1)

include(projectdir("notebooks", "paramsRaw.jl"))

#%%

derivedParams = derivedParameters(Ω, Ωperp, N, k_Cd, k_Ca, k_Sd, k_Sa, k₁, k₂, k₃, k₄, 𝒞, 𝒮, ℰ, D_C, D_S, Tᵣstar; checks=true)
@unpack L₀, E₀, C_b, S_b, δ_C, δ_S, α_C, α_S, C₀, S₀, Tᵣ, T̃ᵣ, K₂, K₃, K₄, σ, ϵ, 𝒟, β, h_C, h_S = derivedParams

#%%

# Create directory for run data labelled with current time.
paramsName = @savename nSpatialDims K₂ K₄ α_C β 𝒟 T̃ᵣ thicknessProfile differencing
folderName = "$(Dates.format(Dates.now(),"yy-mm-dd-HH-MM-SS"))_$(paramsName)"
# Create frames subdirectory to store system state at each output time
mkpath(datadir("sims",subFolder,folderName))

#%%

sol, p = glycosylationAnyD(dims, K₂, K₄, T̃ᵣ, α_C, 𝒟, β, thickness=thicknessProfile, differencing=differencing, solver=solver, nOutputs=nOutputs, σGRF=σGRF, terminateAt=terminateAt)
println("finished sim")

# u = h₀/h_C
# λ = h_C/h_S
# ζ = (2*k₂*Ωperp)/(k₃*𝒮)
# γ = (2*k₂*Ωperp)/(k₁*𝒞)
# Δ = 2*k₂*k₄*Ωperp/(k₁*k₃*𝒮)
# F = (u*(1-Δ*(1+λ*u)))/((1+u)*(1+ζ*(1+λ*u)*(1+u+(1/γ))))

T̃ᵣ₅₀ = sol.t[end]
Tᵣ₅₀ = T̃ᵣ₅₀*(N^2)*(K₂+σ*K₃)
Tᵣ₅₀Star = Tᵣ₅₀/(k₁*E₀)
@show P_star(sol.u[end], p.W, p.dims, p.dν, p.hᵥ, α_C, C_b, Ω, ϕ, Ωperp, k₁, ℰ, Tᵣ₅₀Star)
@show P_analy = Pstar₅₀Analytic(h₀, h_C, h_S, k₁, k₂, k₃, k₄, Ωperp, 𝒮, 𝒞, ℰ, N, ϕ)


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
    𝒞 = 𝒞,
    𝒮 = 𝒮,
    ℰ = ℰ,
    D_C = D_C,
    D_S = D_S,
    Tᵣstar = Tᵣstar,
    ϕ = ϕ
)

jldsave(datadir("sims",subFolder,folderName,"solution.jld2"); sol, p, rawParams)

#%%

if nSpatialDims==1
    concentrationSurfaceMovie(sol.u, dims; subFolder=subFolder, folderName=folderName)
    # concentrationHeatmapMovie(sol.u, dims; subFolder=subFolder, folderName=folderName)
    # spaceIntegralOver_ν_Movie(sol.u, p; subFolder=subFolder, folderName=folderName)
    # if thicknessProfile=="GRF"
    #     thicknessPlot(p.hᵥ, p.dims; subFolder=subFolder, folderName=folderName)
    # end
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



#%%

# h₀ = 0.002
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
# Tᵣstar= 1000000000000.0  # Release time
# ϕ     = 0.5