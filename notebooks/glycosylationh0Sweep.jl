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

thicknessProfile = "uniform"
differencing = "centre"
solver = SSPRK432()
nOutputs = 100

include(projectdir("notebooks", "paramsRaw.jl"))

nSpatialDims = 1
Ngrid = 401
dims = fill(Ngrid, nSpatialDims+1)

h_C = 2*k_Ca/k_Cd
h_S = 2*k_Sa/k_Sd
hMax = h_C*10
hMin = h_C/10
h₀s = collect(hMin:5*hMin:hMax)
Ωs    = h₀s.*Ωperp      # Dimensional lumen volume 


λ = h_C/h_S
ζ = (2*k₂*Ωperp)/(k₃*𝒮)
γ = (2*k₂*Ωperp)/(k₁*𝒞)
Δ = 2*k₂*k₄*Ωperp/(k₁*k₃*𝒮)
F = (u*(1-Δ*(1+λ*u)))/((1+u)*(1+ζ*(1+λ*u)*(1+u+(1/γ))))
hMax = sqrt((γ+ζ)/(ζ*λ))
# 𝒞_bConstant = 𝒞/Ω
# 𝒮_bConstant = 𝒮/Ω
# 𝒞s  = 𝒞_bConstant.*Ωs
# 𝒮s  = 𝒮_bConstant.*Ωs

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

PstarsAnalytic = []
PstarsSim = []
MstarsPhiSim = []
Tᵣ₅₀Stars = []
sols = []
ps = []
for i=1:length(h₀s)
    @show h₀s[i]
    # derivedParams = derivedParameters(Ωs[i], Ωperp, N, k_Cd, k_Ca, k_Sd, k_Sa, k₁, k₂, k₃, k₄, 𝒞s[i], 𝒮s[i], ℰ, D_C, D_S, Tᵣstar; checks=false)
    derivedParams = derivedParameters(Ωs[i], Ωperp, N, k_Cd, k_Ca, k_Sd, k_Sa, k₁, k₂, k₃, k₄, 𝒞, 𝒮, ℰ, D_C, D_S, Tᵣstar; checks=false)
    @unpack L₀, E₀, C_b, S_b, δ_C, δ_S, α_C, α_S, C₀, S₀, Tᵣ, T̃ᵣ, K₂, K₃, K₄, σ, ϵ, 𝒟, β, h_C, h_S = derivedParams
    sol, p = glycosylationAnyD(dims, K₂, K₄, T̃ᵣ, α_C, 𝒟, β, thickness=thicknessProfile, differencing=differencing, solver=solver, nOutputs=nOutputs, terminateAt="halfProduction")
    T̃ᵣ₅₀ = sol.t[end]
    Tᵣ₅₀ = T̃ᵣ₅₀*(N^2)*(K₂+σ*K₃)
    Tᵣ₅₀Star = Tᵣ₅₀/(k₁*E₀)
    push!(Tᵣ₅₀Stars, Tᵣ₅₀Star)
    push!(sols, sol)
    push!(ps, p)
    push!(MstarsPhiSim, M_star_ϕ(sol.u[end], p.W, p.dims, p.dν, p.hᵥ, α_C, 𝒞, Ωs[i], ϕ))
    push!(PstarsSim, P_star(sol.u[end], p.W, p.dims, p.dν, p.hᵥ, α_C, C_b, Ωs[i], ϕ, Ωperp, k₁, ℰ, Tᵣ₅₀Star))
    @show PstarsSim[end]
    push!(PstarsAnalytic, Pstar₅₀Analytic(h₀s[i], h_C, h_S, k₁, k₂, k₃, k₄, Ωperp, 𝒮, 𝒞, ℰ, N, ϕ))
end

#%%

linesVec = []
labelsVec = []
fig = Figure()#size=(500,500))

ax1 = Axis(fig[1,1])
push!(linesVec, lines!(ax1, h₀s, PstarsSim, color=:blue))
push!(labelsVec, "Numerical")
ylims!(ax1, (0.0, maximum(PstarsSim)))
xlims!(ax1, (0.0, maximum(h₀s)))
ax1.xlabel = "h₀"
ax1.ylabel = L"𝓟^*_{50}"

ax2 = Axis(fig[2,1])
ax2.xlabel = "h₀"
ax2.ylabel = L"𝓟^*_{50}"
ylims!(ax2, (0.0,maximum(PstarsAnalytic)))
xlims!(ax2, (0.0, max(maximum(h₀s), max(h_C, h_S))))
push!(linesVec, lines!(ax2, h₀s, PstarsAnalytic, color=:red))
push!(labelsVec, "Analytic")

push!(linesVec, vlines!(ax1, h_C, color=:green))
push!(labelsVec, L"h_C")
push!(linesVec, vlines!(ax1, h_S, color=:orange))
push!(labelsVec, L"h_S")

push!(linesVec, vlines!(ax2, h_C, color=:green))
push!(labelsVec, L"h_C")
push!(linesVec, vlines!(ax2, h_S, color=:orange))
push!(labelsVec, L"h_S")

Legend(fig[:,2], linesVec, labelsVec)

linkxaxes!(ax1, ax2)

display(fig)

subFolder = "h0sweep"
mkpath(datadir("sims", subFolder))
save(datadir("sims", subFolder, "simulationPvsh.png"),fig)

#%%

for i=1:length(sols)
    folderName = "h_0=$(h₀s[i])"
    concentrationSurfaceMovie(sols[i].u, dims; subFolder=subFolder, folderName=folderName) 
end



# Ωperp = 10000    # Dimensional lumen footprint area
# # Ω     = h₀*Ωperp      # Dimensional lumen volume 
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
# 𝓢     = 100000.0
# 𝓔     = 0.0001
# D_C   = 0.0000001  # Monomer/polymer diffusivity
# D_S   = 0.0000001  # Substrate diffusivity
# Tᵣstar= 10000000000000.0  # Release time
# ϕ     = 0.5