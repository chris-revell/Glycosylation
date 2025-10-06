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
nOutputs = 1000

include(projectdir("notebooks", "paramsRaw.jl"))

nSpatialDims = 1
Ngrid = 201
dims = fill(Ngrid, nSpatialDims+1)

h_C = 2*k_Ca/k_Cd
h_S = 2*k_Sa/k_Sd
hMax = h_C*5
hMin = h_C/10
h₀s = collect(hMin:2*hMin:hMax)
Ωs = h₀s.*𝒜      # Dimensional lumen volume 

#%%

PstarsAnalytic = []
PstarsSim = []
PstarsSim50 = []
MstarsPhiSim = []
T̃ᵣ₅₀s = []
Tᵣ₅₀s = []
Tᵣ₅₀Stars = []
sols = []
ps = []
α_Cs = []
for i=1:length(h₀s)
    @show h₀s[i]    
    derivedParams = derivedParameters(Ωs[i], 𝒜, N, k_Cd, k_Ca, k_Sd, k_Sa, k₁, k₂, k₃, k₄, 𝒞, 𝒮, ℰ, D_C, D_S, Tᵣstar; checks=false)
    @unpack L₀, E₀, C_b, S_b, δ_C, δ_S, α_C, α_S, C₀, S₀, Tᵣ, T̃ᵣ, K₂, K₃, K₄, σ, ϵ, 𝒟, β, h_C, h_S, λ, ζ, γ, Δ, F = derivedParams
    
    sol, p = glycosylationAnyD(dims, K₂, K₄, 1000.0, α_C, 𝒟, β, thickness=thicknessProfile, differencing=differencing, solver=solver, nOutputs=nOutputs, terminateAt="halfProduction", saveIntermediate=false)
    
    push!(α_Cs, α_C)
    T̃ᵣ₅₀ = sol.t[end]
    Tᵣ₅₀ = T̃ᵣ₅₀*(N^2)*(K₂+σ*K₃)
    Tᵣ₅₀Star = Tᵣ₅₀/(k₁*E₀)
    push!(sols, sol)
    push!(ps, p)
    push!(T̃ᵣ₅₀s, T̃ᵣ₅₀)
    push!(Tᵣ₅₀s, Tᵣ₅₀)
    push!(Tᵣ₅₀Stars, Tᵣ₅₀Star)
    push!(MstarsPhiSim, Mstarϕ(sol.u[end], p.W, p.dims, p.dν, p.hᵥ, α_C, 𝒞, ϕ))
    push!(PstarsSim50, 𝒫star₅₀Numeric(k₁, k₂, k₃, 𝒞, ℰ, 𝒮, h₀, k_Ca, k_Cd, k_Sa, k_Sd, 𝒜, T̃ᵣ₅₀))
    push!(PstarsSim, Mstarϕ(sol.u[end], p.W, p.dims, p.dν, p.hᵥ, α_C, 𝒞, ϕ)/Tᵣ₅₀Star)
    push!(PstarsAnalytic, 𝒫star₅₀Analytic(h₀s[i], h_C, h_S, k₁, k₂, k₃, k₄, 𝒜, 𝒮, 𝒞, ℰ, N, ϕ))
end

#%%

linesVec = []
labelsVec = []
fig = Figure()#size=(500,500))

ax1 = Axis(fig[1,1])
push!(linesVec, lines!(ax1, h₀s, PstarsSim, color=:blue))
push!(labelsVec, "Numerical")
push!(linesVec, lines!(ax1, h₀s, PstarsAnalytic, color=:red))
push!(labelsVec, "Analytic")
push!(linesVec, vlines!(ax1, h_C, color=:green))
push!(labelsVec, L"h_C")
push!(linesVec, vlines!(ax1, h_S, color=:orange))
push!(labelsVec, L"h_S")

# ylims!(ax1, (0.0, maximum(PstarsSim)))
xlims!(ax1, (0.0, 1.1*max(maximum(h₀s), max(h_C, h_S))))

ax1.xlabel = "h₀"
ax1.ylabel = L"𝓟^*_{50}"

Legend(fig[:,2], linesVec, labelsVec)

display(fig)

# Create directory for run data labelled with current time.
subFolder = "h0sweep"
# paramsName = @savename nSpatialDims K₂ K₄ α_C β 𝒟 T̃ᵣ thicknessProfile differencing
folderName = "$(Dates.format(Dates.now(),"yy-mm-dd-HH-MM-SS"))"#_$(paramsName)"
# Create frames subdirectory to store system state at each output time
mkpath(datadir("sims",subFolder,folderName))

save(datadir("sims", subFolder, folderName, "simulationPvsh.png"), fig)

#%%

for i=1:length(sols)
    # mkpath(datadir("sims",subFolder,"$folderName/h_0=$(h₀s[i])"))
    concentrationSurfaceMovie(sols[i].u, dims; subFolder=datadir("sims",subFolder,folderName), folderName="h_0=$(h₀s[i])") 
end


#%%

fig = Figure(size=(1000, 1000), fontsize=16)
Label(fig[1, 2, Top()], L"Dimensional")
ax1 = Axis(fig[1,2])
ax2 = Axis(fig[2,2])
ax3 = Axis(fig[3,2])
ax4 = Axis(fig[4,2])
ax1.xlabel = L"h_0"
ax2.xlabel = L"h_0"
ax3.xlabel = L"h_0"
ax4.xlabel = L"h_0"
ax1.ylabel = L"M^*_{50}"
ax2.ylabel = L"T^*_{50}"
ax3.ylabel = L"\mathcal{P}^*_{50}"
ax4.ylabel = L"\mathcal{P}^*_{50} (analytic)"
lines!(ax1, h₀s, MstarsPhiSim)
lines!(ax2, h₀s, Tᵣ₅₀Stars)
lines!(ax3, h₀s, PstarsSim)
lines!(ax4, h₀s, PstarsAnalytic)

Label(fig[1, 1, Top()], L"Dimensionless")
ax12 = Axis(fig[1,1])
ax22 = Axis(fig[2,1])
ax32 = Axis(fig[3,1])
# ax42 = Axis(fig[4,1])
ax12.xlabel = L"h_0"
ax22.xlabel = L"h_0"
ax32.xlabel = L"h_0"
# ax42.xlabel = L"h_0"
ax12.ylabel = L"\tilde{M}_{50}"
ax22.ylabel = L"\tilde{T}_{50}"
ax32.ylabel = L"\mathcal{\tilde{P}}_{50}"
# ax42.ylabel = L"\mathcal{P}^*_{50} (analytic)"
lines!(ax12, h₀s, fill(π/2.0, length(h₀s)))
lines!(ax22, h₀s, T̃ᵣ₅₀s)
lines!(ax32, h₀s, fill(π/2.0, length(h₀s))./T̃ᵣ₅₀s)
# lines!(ax4, h₀s, PstarsAnalytic)
display(fig)
save(datadir("sims", subFolder, folderName, "fullbreakdown.png"), fig)
