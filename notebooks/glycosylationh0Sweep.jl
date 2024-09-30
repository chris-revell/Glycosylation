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
using CairoMakie

@from "$(srcdir("Glycosylation.jl"))" using Glycosylation
@from "$(srcdir("Visualise.jl"))" using Visualise
@from "$(srcdir("UsefulFunctions.jl"))" using UsefulFunctions
@from "$(srcdir("MakeWeightMatrices.jl"))" using MakeWeightMatrices
@from "$(srcdir("DerivedParameters.jl"))" using DerivedParameters

nSpatialDims = 1

Ωperp = 1.0  # Lumen footprint area
N     = 100         # Maximum polymer length 
k_Cd  = 3000.0 # Complex desorption rate
k_Ca  = 0.01 # Complex adsorption rate
k_Sd  = 1.0 # Substrate desorption rate
k_Sa  = 1.0 # Substrate adsorption rate
k₁    = 2.0   # Complex formation forward reaction rate 
k₂    = 0.01   # Complex dissociation reverse reaction rate 
k₃    = 0.01   # Product formation
k₄    = 2.0  # Product dissociation 
E_0   = 0.01
𝓒     = 1.0
𝓢     = 1000.0
D_C   = 0.000001  # Monomer/polymer diffusivity
D_S   = 0.000001  # Substrate diffusivity
Tᵣstar= 0.1  # Release time
ϕ     = 0.5

Ngrid = 101
nSpatialDims == 1 ? dims  = [Ngrid, Ngrid] : dims  = [Ngrid, Ngrid, Ngrid]
xMax = (Ωperp)^(1/nSpatialDims)
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

W = vertexVolumeWeightsMatrix(dims, spacing)

# h₀s = collect(0.1:0.1:3.0)
h₀s = collect(0.000001:0.000001:0.00010)

sols = []
hᵥs = []
α_Cs = []
Ωs =[]
C_bs =[]

for h₀ in h₀s
    @show h₀
    mat_h = h₀.*ones(fill(Ngrid, nSpatialDims+1)...)
    derivedParams = derivedParameters(h₀, Ωperp, N, k_Cd, k_Ca, k_Sd, k_Sa, k₁, k₂, k₃, k₄, E_0, 𝓒, 𝓢, D_C, D_S, Tᵣstar; checks=false)
    @unpack K₂, K₄, Tᵣ, α_C, 𝓓, β, Ω, C_b = derivedParams
    sol = glycosylationAnyD(mat_h, dims, Ωperp, 𝓒, K₂, K₄, Tᵣ, α_C, 𝓓, β)
    push!(sols, sol)
    push!(Ωs, Ω)
    push!(α_Cs, α_C)
    push!(C_bs, C_b)
    hᵥ_vec = reshape(mat_h, (Ngrid)^(nSpatialDims+1))
    hᵥ = spdiagm(hᵥ_vec)
    push!(hᵥs, hᵥ)
end

#%%


h_C = 2*k_Ca/k_Cd
h_S = 2*k_Sa/k_Sd


fig = Figure(size=(500,500))
ax1 = Axis(fig[1,1])
PstarsSim = Float64[]
for i=1:length(sols)
    push!(PstarsSim, P_star(sols[i][end], W, dims, hᵥs[i], ϕ, α_Cs[i], C_bs[i], Ωs[i], spacing[1], Tᵣstar))
end
linesVec = []
labelsVec = []
push!(linesVec, lines!(ax1, h₀s, PstarsSim, color=:blue))
push!(labelsVec, "Numerical")
ylims!(ax1, (0.0,maximum(PstarsSim)))
# xlims!(ax1, (0.0,maximum(h₀s)))
ax1.xlabel = "h₀"
ax1.ylabel = L"𝓟^*"

PstarsAnalytic = 𝓟starUniform.(ϕ, 𝓒, 𝓢, E_0, h₀s, Ωperp, k_Ca, k_Cd, k_Sa, k_Sd, k₁, k₂, k₃, k₄, N, Tᵣstar)

ax2 = Axis(fig[2,1])
ax2.xlabel = "h₀"
ax2.ylabel = L"𝓟^*"
ylims!(ax2, (0.0,maximum(PstarsAnalytic)))
# xlims!(ax2, (0.0,maximum(h₀s)))
push!(linesVec, lines!(ax2, h₀s, PstarsAnalytic, color=:red))
push!(labelsVec, "Analytic")
# lines!(ax1, h₀s, PstarsAnalytic, label="Analytic")



# linkxaxes!(ax1, ax2)

push!(linesVec, vlines!(ax1, h_C, color=:green))
push!(labelsVec, L"h_C")
push!(linesVec, vlines!(ax2, h_C, color=:green))
push!(labelsVec, L"h_C")
# push!(linesVec, vlines!(ax1, h_S, color=:orange))
# push!(labelsVec, L"h_S")
# push!(linesVec, vlines!(ax2, h_S, color=:orange))
# push!(labelsVec, L"h_S")

Legend(fig[:,2], linesVec, labelsVec)

display(fig)

# save("simulationPvsh.png",fig)



