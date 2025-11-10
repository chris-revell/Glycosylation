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
using Statistics
using JLD2

@from "$(srcdir("Glycosylation.jl"))" using Glycosylation
@from "$(srcdir("UsefulFunctions.jl"))" using UsefulFunctions
@from "$(srcdir("MakeWeightMatrices.jl"))" using MakeWeightMatrices
@from "$(srcdir("DerivedParameters.jl"))" using DerivedParameters

#%%

subFolder = "Figure3"
folderName = "25-10-30-17-26-58"
data1 = load(datadir("sims", subFolder, folderName, "solutions.jld2"))
@unpack sols, ps, h₀s = data1

include(projectdir("notebooks", "paramsRaw.jl"))
h_C = 2*k_Ca/k_Cd
h_S = 2*k_Sa/k_Sd
Ωs = h₀s.*𝒜

#%%

𝒫sim = []
𝒫analytic = [0.0]
𝒫analyticAdjusted = [0.0]
for i=1:length(h₀s)
    derivedParams = derivedParameters(Ωs[i], 𝒜, N, k_Cd, k_Ca, k_Sd, k_Sa, k₁, k₂, k₃, k₄, 𝒞, 𝒮, ℰ, D_C, D_S, Tᵣstar; checks=false)
    @unpack L₀, E₀, C_b, S_b, δ_C, δ_S, α_C, α_S, C₀, S₀, Tᵣ, T̃ᵣ, K₂, K₃, K₄, σ, ϵ, 𝒟, β, h_C, h_S, λ, ζ, γ, Δ, F = derivedParams
    Tᵣ₅₀Star = sols[i].t[end]*(N^2)*(K₂+σ*K₃)/(k₁*E₀)
    Tᵣ₅₀Analytic = ((N^2)*(K₂+σ*K₃)/(k₁*E₀))*T̃ᵣ₅₀Analytic(𝒞, 𝒮, ϕ, h₀s[i], h_C, h_S, k₁, k₂, k₃, k₄, 𝒜, N, 0.0, 0.0)
    push!(𝒫sim, Mstarϕ(sols[i].u[end], ps[i].W, ps[i].dims, ps[i].dν, ps[i].hᵥ, α_C, 𝒞, ϕ)/Tᵣ₅₀Star)
    push!(𝒫analytic, 𝒫star₅₀Analytic(h₀s[i], h_C, h_S, k₁, k₂, k₃, k₄, 𝒜, 𝒮, 𝒞, ℰ, N, ϕ) )

    midpoint = length(sols[i].u)÷2
    C_peak, ind_peak = findmax(reshape(sols[i].u[midpoint], ps[i].dims...)[:,1])
    νs   = collect(range(0.0, 1.0, ps[i].dims[1]))
    ν_peak = νs[ind_peak]
    Ẽ = ps[i].K₂/(1+ps[i].K₂)
    D = Ẽ*ps[i].K₂*K₄/(1+α_C)
    t̃₀ = sols[i].t[midpoint] - 1/(4.0*π*D*C_peak^2)
    ν₀ = ν_peak - Ẽ*β*(sols[i].t[midpoint]-t̃₀)/(1+α_C)
    Tᵣ₅₀AnalyticAdjusted = ((N^2)*(K₂+σ*K₃)/(k₁*E₀))*T̃ᵣ₅₀Analytic(𝒞, 𝒮, ϕ, h₀s[i], h_C, h_S, k₁, k₂, k₃, k₄, 𝒜, N, ν₀, t̃₀)
    push!(𝒫analyticAdjusted, (α_C*𝒞/(π*(1+α_C)))*(π/2)/Tᵣ₅₀AnalyticAdjusted )
end

ς₅₀ = []
h₀s2 = collect(range(minimum(h₀s), maximum(h₀s), 100))
for h₀ in h₀s2
    derivedParams = derivedParameters(𝒜*h₀, 𝒜, N, k_Cd, k_Ca, k_Sd, k_Sa, k₁, k₂, k₃, k₄, 𝒞, 𝒮, ℰ, D_C, D_S, Tᵣstar; checks=false)
    @unpack L₀, E₀, C_b, S_b, δ_C, δ_S, α_C, α_S, C₀, S₀, Tᵣ, T̃ᵣ, K₂, K₃, K₄, σ, ϵ, 𝒟, β, h_C, h_S, λ, ζ, γ, Δ, F = derivedParams
    push!(ς₅₀, 2*K₂*K₄*ϕ/β)
end

hcutoff = (2.0*k_Sa/k_Sd)*((𝒮*k₁*k₃)/(2.0*𝒜*k₂*k₄) - 1.0)
for h₀cut in [hcutoff]
    derivedParams = derivedParameters(𝒜*h₀cut, 𝒜, N, k_Cd, k_Ca, k_Sd, k_Sa, k₁, k₂, k₃, k₄, 𝒞, 𝒮, ℰ, D_C, D_S, Tᵣstar; checks=false)
    @unpack L₀, E₀, C_b, S_b, δ_C, δ_S, α_C, α_S, C₀, S₀, Tᵣ, T̃ᵣ, K₂, K₃, K₄, σ, ϵ, 𝒟, β, h_C, h_S, λ, ζ, γ, Δ, F = derivedParams

    # midpoint = length(sol1.u)
    # C_peak, ind_peak = findmax(reshape(sol1.u[midpoint], p1.dims...)[:,1])
    # νs   = collect(range(0.0, 1.0, p1.dims[1]))
    # ν_peak = νs[ind_peak]
    # Ẽ = p1.K₂/(1+p1.K₂)
    # D = Ẽ*p1.K₂*K₄/(1+α_C)
    # t̃₀ = sol1.t[midpoint] - 1/(4.0*π*D*C_peak^2)
    # ν₀ = ν_peak - Ẽ*β*(sol1.t[midpoint]-t̃₀)/(1+α_C)
    
    Tᵣ₅₀Analytic = ((N^2)*(K₂+σ*K₃)/(k₁*E₀))*T̃ᵣ₅₀Analytic(𝒞, 𝒮, ϕ, h₀cut, h_C, h_S, k₁, k₂, k₃, k₄, 𝒜, N, 0.0, 0.0)
    # push!(𝒫analytic, (α_C*𝒞/(π*(1+α_C)))*(π/2)/Tᵣ₅₀Analytic )
    push!(𝒫analytic, 𝒫star₅₀Analytic(h₀cut, h_C, h_S, k₁, k₂, k₃, k₄, 𝒜, 𝒮, 𝒞, ℰ, N, ϕ) )
    push!(𝒫analyticAdjusted, 0.0 )
end

#%%

differencing = "centre"
solver = SSPRK432()
nOutputs = 1000
nSpatialDims = 1
Ngrid = 401
dims = fill(Ngrid, nSpatialDims+1)

include(projectdir("notebooks", "paramsRaw.jl"))

derivedParams = derivedParameters(Ω, 𝒜, N, k_Cd, k_Ca, k_Sd, k_Sa, k₁, k₂, k₃, k₄, 𝒞, 𝒮, ℰ, D_C, D_S, Tᵣstar; checks=false)
@unpack L₀, E₀, C_b, S_b, δ_C, δ_S, α_C, α_S, C₀, S₀, Tᵣ, T̃ᵣ, K₂, K₃, K₄, σ, ϵ, 𝒟, β, h_C, h_S, λ, ζ, γ, Δ, F = derivedParams

sol1, p1 = glycosylation(dims, K₂, K₄, T̃ᵣ, α_C, 𝒟, β, thickness="uniform", differencing="centre", solver=solver, nOutputs=nOutputs, terminateAt="nuWall")
println("finished sim 1")
sol2, p2 = glycosylation(dims, K₂, K₄, T̃ᵣ, α_C, 𝒟, β, thickness="uniform", differencing=differencing, solver=solver, nOutputs=nOutputs, terminateAt="halfProduction", saveIntermediate=false) 
println("finished sim 2")
Tᵣ₅₀Star = sol2.t[end]*(N^2)*(K₂+σ*K₃)/(k₁*E₀)
𝒫sim1 = Mstarϕ(sol2.u[end], p2.W, p2.dims, p2.dν, p2.hᵥ, α_C, 𝒞, 0.5)/Tᵣ₅₀Star
@show 𝒫sim1

#%%

stoppoint = -5 #-15

midpoint = (length(sol1.u)+stoppoint)÷2
C_peak, ind_peak = findmax(reshape(sol1.u[midpoint], p1.dims...)[:,p1.dims[2]÷2])
νs   = collect(range(0.0, 1.0, p1.dims[1]))
ν_peak = νs[ind_peak]
Ẽ = p1.K₂/(1+p1.K₂)
D = Ẽ*p1.K₂*K₄/(1+α_C)
t̃₀ = sol1.t[midpoint] - 1/(4.0*π*D*C_peak^2)
ν₀ = ν_peak - Ẽ*β*(sol1.t[midpoint]-t̃₀)/(1+α_C)
νsOffset = νs.-ν₀
tsOffset = sol1.t.-t̃₀
firstPositivetIndex = findfirst(x->x>0, tsOffset)

fig = Figure(size=(1000,800), fontsize=24, figure_padding=(0,50,1,11))
g1 = GridLayout(fig[1,1])
g2 = GridLayout(fig[2,1])

ax1 = Axis(g1[1, 1])
allLines_ax1 = []
allTs_ax1 = []
colorsUsed = [(:red), (:green), (:blue)]
for (c,i) in enumerate([firstPositivetIndex, (length(sol1.t)+stoppoint-firstPositivetIndex)÷2+firstPositivetIndex, length(sol1.t)+stoppoint])
    uInternal = reshape(sol1.u[i], p1.dims...)
    push!(allLines_ax1, lines!(ax1, νs, uInternal[:,1], linestyle=:solid, color=(colorsUsed[c], 0.5), linewidth=4))
    push!(allLines_ax1, lines!(ax1, νs, homogeneousWidthC.(νsOffset, p1.K₂, K₄, α_C, β, tsOffset[i]), linestyle=:dot, color=(colorsUsed[c], 1.0), linewidth=4))
    push!(allTs_ax1, @sprintf("%.2f", tsOffset[i]))
end
labels_ax1 = []
for t in allTs_ax1
    push!(labels_ax1, "Numeric, t=$t")
    push!(labels_ax1, "Asymptotic, t=$t")
end
ax1.xlabel = L"\nu"
ax1.ylabel = L"\tilde{C}"
ax1.yticks = (0.0:20.0:20.0, [L"0.0", L"20.0"])
ax1.xticks = (0.0:0.5:1.0, [L"0.0", L"\nu=\phi", L"1.0"])
lϕ = vlines!(ax1, 0.5, color=(:black, 0.5))#, linewidth=4
text!(ax1, Point{2,Float64}(0.95,0.9*20.0), text="A", color=:black, align=(:right, :bottom), fontsize=24) 
text!(ax1, Point{2,Float64}(0.05, 0.9*20.0), text = L"\tilde{t}=%$(allTs_ax1[1])", color=:red) 
text!(ax1, Point{2,Float64}(0.23, 0.38*20.0), text = L"\tilde{t}=%$(allTs_ax1[2])", color=:green) 
text!(ax1, Point{2,Float64}(0.65, 0.28*20.0), text = L"\tilde{t}=%$(allTs_ax1[3])", color=:blue) 
ylims!(ax1, (0.0, 20.0))
xlims!(ax1, (0.0, 1.0))


ax2 = Axis(g1[1,2], yticks = (0.0:π/2.0:π, [L"0", L"π/2", L"π"]))#, aspect=AxisAspect(1.5))
ylims!(ax2, (0.0, 1.05*π))
xlims!(ax2, (0.0, sol1.t[end+stoppoint]))
tSeries = sol1.t[firstPositivetIndex:end+stoppoint]
numericalMs = [M̃ϕ(u, p1.W, p1.dims, p1.dν, p1.hᵥ, ϕ, thresh="ceil") for u in sol1.u[firstPositivetIndex:end+stoppoint]]
analyticMs = [M̃ϕAnalytic.(ϕ, ν₀, τ, α_C, β, p1.K₂, K₄) for τ in tsOffset[firstPositivetIndex:end+stoppoint]]
allLines = [ lines!(ax2, tSeries, numericalMs, linewidth=4, color=(:red, 0.5)), 
                lines!(ax2, tSeries, analyticMs, linewidth=4, color=(:blue, 1.0), linestyle=:dot),                 
            ]
allLabels = [ "Numeric",
    "Asymptotic",
]
ind = findfirst(x->M̃ϕ(x, p1.W, p1.dims, p1.dν, p1.hᵥ, ϕ)>=π/2.0, sol1.u)
l50 = vlines!(ax2, sol1.t[ind], color=(:black, 0.5))#, linewidth=4
tEndString = @sprintf("%.2f", sol1.t[end+stoppoint])
ax2.xticks = ([0.0, sol1.t[ind], sol1.t[end+stoppoint]], [L"0.0", L"\tilde{t}=\tilde{T}_{r,50}", L"%$(tEndString)"])
ax2.xlabel = L"\tilde{t}"
ax2.ylabel = L"\tilde{M}_\phi"
text!(ax2, Point{2,Float64}(0.95*sol1.t[end+stoppoint],0.9*π), text="B", color=:black, align=(:right, :bottom), fontsize=24) 


ax3 = Axis(g2[1,1])
linesVec_ax3 = []
labelsVec_ax3 = []
hcutoff = (2.0*k_Sa/k_Sd)*((𝒮*k₁*k₃)/(2.0*𝒜*k₂*k₄) - 1.0)
push!(linesVec_ax3, lines!(ax3, h₀s, 𝒫sim, color=(:red, 0.5), linewidth=4))
push!(labelsVec_ax3, "Numeric")
push!(linesVec_ax3, lines!(ax3, [0.0, h₀s..., hcutoff], 𝒫analytic, color=(:blue, 1.0), linewidth=4, linestyle=:dot))
push!(labelsVec_ax3, "Asymptotic")
push!(linesVec_ax3, vlines!(ax3, h_C, color=(:black, 0.5)))#, linewidth=4))
push!(linesVec_ax3, vlines!(ax3, h_S, color=(:black, 0.5)))#, linewidth=4))
push!(linesVec_ax3, vlines!(ax3, hcutoff, color=(:black, 0.5)))#, linewidth=4))
# ax3.xticks = ([0.0, h_C, h_S, hcutoff], [L"0", L"  h_C", L"h_S", L"h_{cut-off}"])
ax3.xticks = ([0.0, h_C, h_S, 3.0, 6.0, hcutoff], [L"0", L"h_C", L"h_S", L"3.0", L"6.0", L"h_{cut-off}"])
# ax3.xticks = ([0.0, h_C, 1.0, h_S], [L"0", L"h_C", L"1.0", L"h_S"])
ax3.yticks = ([0.0, 0.0001], [L"0.0", L"0.0001"])
ax3.xaxis.elements[:ticklabels].align = tuple.([:right, :left, :center, :center, :center, :center], :top)
ax3.xgridvisible=false
xlims!(ax3, (0.0, 1.05*hcutoff))
ylims!(ax3, (0.0, 1.1*maximum(𝒫analytic)))
ax3.xlabel = L"h_0 [L]"
ax3.ylabel = L"𝓟^*_{50} [M/T]"
text!(ax3, Point{2,Float64}(0.9*1.05*hcutoff, 0.9*1.1*maximum(𝒫analytic)), text="C", color=:black, align=(:right, :bottom), fontsize=24) 
scatter!(ax3, [1.0], [𝒫sim1], color=:black, marker=:cross, markersize=15)

ax4 = Axis(g2[1,1], yticklabelcolor = :green, ylabelcolor=:green, yaxisposition = :right)
hidespines!(ax4)
hidexdecorations!(ax4)
linesVec_ax4 = []
labelsVec_ax4 = []
hcutoff = (2.0*k_Sa/k_Sd)*((𝒮*k₁*k₃)/(2.0*𝒜*k₂*k₄) - 1.0)
push!(linesVec_ax4, lines!(ax4, h₀s2, sqrt.(ς₅₀), color=(:green, 0.5), linewidth=4))
# push!(labelsVec_ax4, "Numeric")
# push!(linesVec_ax4, lines!(ax4, [0.0, h₀s..., hcutoff], 𝒫analytic, color=(:blue, 1.0), linewidth=4, linestyle=:dot))
# push!(labelsVec_ax4, "Asymptotic")
push!(linesVec_ax4, vlines!(ax4, h_C, color=(:black, 0.5)))#, linewidth=4))
push!(linesVec_ax4, vlines!(ax4, h_S, color=(:black, 0.5)))#, linewidth=4))
push!(linesVec_ax4, vlines!(ax4, hcutoff, color=(:black, 0.5)))#, linewidth=4))
# ax4.xticks = ([0.0, h_C, h_S, hcutoff], [L"0", L"  h_C", L"h_S", L"h_{cut-off}"])
ax4.xticks = ([0.0, h_C, h_S, 3.0, 6.0, hcutoff], [L"0", L"h_C", L"h_S", L"3.0", L"6.0", L"h_{cut-off}"])
# ax4.xticks = ([0.0, h_C, 1.0, h_S], [L"0", L"h_C", L"1.0", L"h_S"])
# ax4.yticks = ([0.0, 0.0001], [L"0.0", L"0.0001"])
ax4.xaxis.elements[:ticklabels].align = tuple.([:right, :left, :center, :center, :center, :center], :top)
ax4.xgridvisible=false
ax4.ygridvisible=false
xlims!(ax4, (0.0, 1.05*hcutoff))
# ylims!(ax4, (0.0, 0.55))
# ax4.xlabel = L"h_0 [L]"
ax4.ylabel = L"\varsigma _{50} [L]"
# text!(ax4, Point{2,Float64}(0.9*1.05*hcutoff, 0.9*1.1*maximum(𝒫analytic)), text="C", color=:black, align=(:right, :bottom), fontsize=24) 

colsize!(g2, 1, Aspect(1, 1.5))


# @show t̃₀
# @show ν₀


# maxProductionValSim, maxProductionIndSim = findmax(𝒫sim) 
# @show maxProductionValSim
# @show h₀s[maxProductionIndSim]

# maxProductionValAnalytic, maxProductionIndAnalytic = findmax(𝒫analytic) 
# @show maxProductionValAnalytic
# @show h₀s[maxProductionIndAnalytic]

resize_to_layout!(fig)

# save(datadir("sims", subFolder, folderName, "Figure3_ν0=$(@sprintf("%.6f", ν₀))_t0=$(@sprintf("%.6f", t̃₀)).png"), fig)
# save(datadir("sims", subFolder, folderName, "Figure3_ν0=$(@sprintf("%.6f", ν₀))_t0=$(@sprintf("%.6f", t̃₀)).pdf"), fig)
save(datadir("sims", subFolder, folderName, "Figure3.pdf"), fig)
display(fig)