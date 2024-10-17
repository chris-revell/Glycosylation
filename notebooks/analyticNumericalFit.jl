#%%
# flux_νₑ = (diffusive_flux_ν + advective_flux_ν)
# flux_νₑ = K₂*K₄.*Pν*∇ₑ*cᵥ - β*Pν*cₑ    where cᵥ is concentration over vertices, cₑ is concentration over edges 
# cₑ = Aᵤₚ*cᵥ
# flux_νₑ = (K₂*K₄.*Pν*∇ₑ - β*Pν*Aᵤₚ)*cᵥ
# flux_xyₑ = Dₑ*hₑ*diffusive_flux_xy
# flux_xyₑ = Dₑ*hₑ*Pxy*∇ₑ*cᵥ
# ċ = aE∇⋅flux_νₑ + a∇⋅flux_xyₑ
# ċ = a*E*∇⋅(K₂*K₄.*Pν*∇ₑ*cᵥ - β*Pν*Aᵤₚ*cᵥ) + a∇⋅(Dₑ*hₑ*Pxy*∇ₑ*cᵥ)
# Dₑ constant over edges 
# ċ = a*(E*∇⋅(K₂*K₄.*Pν*∇ₑ - β*Pν*Aᵤₚ) + 𝓓.*∇⋅(hₑ*Pxy*∇ₑ))*cᵥ

# L = -W⁻¹*Aᵀ*𝓓*l⁻¹*A .+ W⁻¹*Aᵀ*V*Aᵤₚ # Express model as a matrix operator 


# Cνν = W⁻¹*Aᵀ*Pν*l⁻¹*A
# Cν = Aᵀ*l⁻¹*Pν*Aᵤₚ
# flux_νₑ = (diffusive_flux_ν + advective_flux_ν)
# flux_νₑ = K₂*K₄.*Pν*∇ₑ*cᵥ - β*Pν*cₑ    where cᵥ is concentration over vertices, cₑ is concentration over edges 
# cₑ = Aᵤₚ*cᵥ
# flux_νₑ = (K₂*K₄.*Pν*∇ₑ - β*Pν*Aᵤₚ)*cᵥ
# flux_xyₑ = Dₑ*hₑ*diffusive_flux_xy
# flux_xyₑ = Dₑ*hₑ*Pxy*∇ₑ*cᵥ
# ċ = aE∇⋅flux_νₑ + a∇⋅flux_xyₑ
# ċ = a*E*∇⋅(K₂*K₄.*Pν*∇ₑ*cᵥ - β*Pν*Aᵤₚ*cᵥ) + a∇⋅(Dₑ*hₑ*Pxy*∇ₑ*cᵥ)
# Dₑ constant over edges 
# ċ = a*(E*∇⋅(K₂*K₄.*Pν*∇ₑ - β*Pν*Aᵤₚ) + 𝓓.*∇⋅(hₑ*Pxy*∇ₑ))*cᵥ


#


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


h₀ = 1.0

nSpatialDims = 1

K₂ = 1.0
K₄ = 0.00001
Tᵣ = 3.0
α_C = 1.0
𝓓 = 1.0
β = 1.0

Ngrid = 51
dims  = [Ngrid, Ngrid]
#%%

xMax = sqrt(π)
xs   = collect(range(0.0, xMax, dims[2]))
dx   = xs[2]-xs[1]
νMax = 1.0
νs   = collect(range(0.0, νMax, dims[1])) # Positions of discretised vertices in polymerisation space 
dν   = νs[2]-νs[1]
spacing  = [dν, dx]

mat_h = h₀.*ones(fill(Ngrid, nSpatialDims+1)...)

sol = glycosylationAnyD(mat_h, dims, K₂, K₄, Tᵣ, α_C , 𝓓, β)
println("finished sim")

#%%

# Create directory for run data labelled with current time.
paramsName = @savename nSpatialDims K₂ K₄ α_C β 𝓓 Tᵣ h₀
folderName = "$(Dates.format(Dates.now(),"yy-mm-dd-HH-MM-SS"))_$(paramsName)"
# Create frames subdirectory to store system state at each output time
subFolder = ""
mkpath(datadir("sims",subFolder,folderName))

midpoint = length(sol.u)÷3

C_peak, ind_peak = findmax(reshape(sol.u[midpoint], dims...)[:,15])
ν_peak = νs[ind_peak]

Ẽ = K₂/(1+K₂)
D = Ẽ*K₂*K₄/(1+α_C)

t₀ = sol.t[midpoint] - 1/(4.0*π*D*C_peak^2)

ν₀ = ν_peak - Ẽ*β*(sol.t[midpoint]-t₀)/(1+α_C)


#%%

fig = Figure(size=(1000,1000))
ax = CairoMakie.Axis(fig[1, 1], aspect=1)
ax.xlabel = "ν"
ax.ylabel = "C"
analyticLine = Observable(zeros(dims[1]))
numericLine = Observable(zeros(dims[1]))
l1 = lines!(ax, νs, analyticLine, color=:red)
l2 = lines!(ax, νs, numericLine, color=:blue)
Legend(fig[1,2], [l1, l2], ["Analytic", "Numeric"])
ylims!(ax, (0.0, 15.0))
νsOffset = νs.-ν₀
tsOffset = sol.t.-t₀
analyticVals = homogeneousWidthC.(νsOffset, K₂, K₄, α_C, β, tsOffset[1])
record(fig, datadir("sims",subFolder, folderName, "analyticCs.mp4"), 1:length(sol.t); framerate=50) do i
    analyticVals .= homogeneousWidthC.(νsOffset, K₂, K₄, α_C, β, tsOffset[i])
    analyticLine[] .= tst
    uInternal = reshape(sol.u[i], dims...)
    numericLine[] .= uInternal[:,dims[2]÷2]
    analyticLine[] = analyticLine[]
    numericLine[] = numericLine[]
    if i in [1, 251, 500]
        save(datadir("sims",subFolder, folderName, "analyticCs$i.png"), fig)
    end
end


# tst = homogeneousWidthC.(νs, K₂, K₄, α_C, β, sol.t[5])
# analyticLine[] .= tst
# uInternal = reshape(sol.u[end], dims...)
# numericLine[] .= uInternal[:,dims[2]÷2]
# analyticLine[] = analyticLine[]
# numericLine[] = numericLine[]

# display(fig)