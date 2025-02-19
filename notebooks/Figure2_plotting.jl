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

@from "$(srcdir("Glycosylation.jl"))" using Glycosylation
# @from "$(srcdir("Visualise.jl"))" using Visualise
@from "$(srcdir("UsefulFunctions.jl"))" using UsefulFunctions
@from "$(srcdir("MakeWeightMatrices.jl"))" using MakeWeightMatrices
@from "$(srcdir("DerivedParameters.jl"))" using DerivedParameters

#%%

subFolder = "Figure2"
folderName = "25-02-13-16-34-25"
data1 = load(datadir("sims", subFolder, folderName, "solutions.jld2"))
@unpack sols, ps, h₀s = data1

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

h_C = 2*k_Ca/k_Cd
h_S = 2*k_Sa/k_Sd
hMax = h_C*5
hMin = h_C/10
Ωs = h₀s.*Ωperp

#%%

𝒫sim = []
𝒫analytic = [0.0]
𝒫analyticAdjusted = [0.0]
for i=1:length(h₀s)
    @show h₀s[i]    
    derivedParams = derivedParameters(Ωs[i], Ωperp, N, k_Cd, k_Ca, k_Sd, k_Sa, k₁, k₂, k₃, k₄, 𝒞, 𝒮, ℰ, D_C, D_S, Tᵣstar; checks=false)
    @unpack L₀, E₀, C_b, S_b, δ_C, δ_S, α_C, α_S, C₀, S₀, Tᵣ, T̃ᵣ, K₂, K₃, K₄, σ, ϵ, 𝒟, β, h_C, h_S, λ, ζ, γ, Δ, F = derivedParams
    # sol, p = glycosylationAnyD(dims, K₂, K₄, 1000.0, α_C, 𝒟, β, thickness=thicknessProfile, differencing=differencing, solver=solver, nOutputs=nOutputs, terminateAt="halfProduction", saveIntermediate=false) 
    Tᵣ₅₀Star = sols[i].t[end]*(N^2)*(K₂+σ*K₃)/(k₁*E₀)
    Tᵣ₅₀Analytic = ((N^2)*(K₂+σ*K₃)/(k₁*E₀))*T̃ᵣ₅₀Analytic(𝒞, 𝒮, ϕ, h₀s[i], h_C, h_S, k₁, k₂, k₃, k₄, Ωperp, N, 0.0, 0.0)
    push!(𝒫sim, Mstarϕ(sols[i].u[end], ps[i].W, ps[i].dims, ps[i].dν, ps[i].hᵥ, α_C, 𝒞, ϕ)/Tᵣ₅₀Star)
    # push!(𝒫analytic, (α_C*𝒞/(π*(1+α_C)))*(π/2)/Tᵣ₅₀Analytic )
    push!(𝒫analytic, 𝒫star₅₀Analytic(h₀s[i], h_C, h_S, k₁, k₂, k₃, k₄, Ωperp, 𝒮, 𝒞, ℰ, N, ϕ) )

    midpoint = length(sols[i].u)
    C_peak, ind_peak = findmax(reshape(sols[i].u[midpoint], ps[i].dims...)[:,1])
    νs   = collect(range(0.0, 1.0, ps[i].dims[1]))
    ν_peak = νs[ind_peak]
    Ẽ = ps[i].K₂/(1+ps[i].K₂)
    D = Ẽ*ps[i].K₂*K₄/(1+α_C)
    t̃₀ = sols[i].t[midpoint] - 1/(4.0*π*D*C_peak^2)
    ν₀ = ν_peak - Ẽ*β*(sols[i].t[midpoint]-t̃₀)/(1+α_C)
    Tᵣ₅₀AnalyticAdjusted = ((N^2)*(K₂+σ*K₃)/(k₁*E₀))*T̃ᵣ₅₀Analytic(𝒞, 𝒮, ϕ, h₀s[i], h_C, h_S, k₁, k₂, k₃, k₄, Ωperp, N, ν₀, t̃₀)
    push!(𝒫analyticAdjusted, (α_C*𝒞/(π*(1+α_C)))*(π/2)/Tᵣ₅₀AnalyticAdjusted )
end

hcutoff = (2.0*k_Sa/k_Sd)*((𝒮*k₁*k₃)/(2.0*Ωperp*k₂*k₄) - 1.0)
for h₀cut in [hcutoff]
    derivedParams = derivedParameters(Ωperp*h₀cut, Ωperp, N, k_Cd, k_Ca, k_Sd, k_Sa, k₁, k₂, k₃, k₄, 𝒞, 𝒮, ℰ, D_C, D_S, Tᵣstar; checks=false)
    @unpack L₀, E₀, C_b, S_b, δ_C, δ_S, α_C, α_S, C₀, S₀, Tᵣ, T̃ᵣ, K₂, K₃, K₄, σ, ϵ, 𝒟, β, h_C, h_S, λ, ζ, γ, Δ, F = derivedParams

    # midpoint = length(sol1.u)
    # C_peak, ind_peak = findmax(reshape(sol1.u[midpoint], p1.dims...)[:,1])
    # νs   = collect(range(0.0, 1.0, p1.dims[1]))
    # ν_peak = νs[ind_peak]
    # Ẽ = p1.K₂/(1+p1.K₂)
    # D = Ẽ*p1.K₂*K₄/(1+α_C)
    # t̃₀ = sol1.t[midpoint] - 1/(4.0*π*D*C_peak^2)
    # ν₀ = ν_peak - Ẽ*β*(sol1.t[midpoint]-t̃₀)/(1+α_C)
    
    Tᵣ₅₀Analytic = ((N^2)*(K₂+σ*K₃)/(k₁*E₀))*T̃ᵣ₅₀Analytic(𝒞, 𝒮, ϕ, h₀cut, h_C, h_S, k₁, k₂, k₃, k₄, Ωperp, N, 0.0, 0.0)
    # push!(𝒫analytic, (α_C*𝒞/(π*(1+α_C)))*(π/2)/Tᵣ₅₀Analytic )
    push!(𝒫analytic, 𝒫star₅₀Analytic(h₀cut, h_C, h_S, k₁, k₂, k₃, k₄, Ωperp, 𝒮, 𝒞, ℰ, N, ϕ) )
    push!(𝒫analyticAdjusted, 0.0 )
end

#%%

Ωperp = 10000
h₀ = 1.0
Ω = Ωperp*h₀
derivedParams = derivedParameters(Ω, Ωperp, N, k_Cd, k_Ca, k_Sd, k_Sa, k₁, k₂, k₃, k₄, 𝒞, 𝒮, ℰ, D_C, D_S, Tᵣstar; checks=false)
@unpack L₀, E₀, C_b, S_b, δ_C, δ_S, α_C, α_S, C₀, S₀, Tᵣ, T̃ᵣ, K₂, K₃, K₄, σ, ϵ, 𝒟, β, h_C, h_S, λ, ζ, γ, Δ, F = derivedParams

sol1, p1 = glycosylationAnyD(dims, K₂, K₄, T̃ᵣ*3, α_C, 𝒟, β, thickness=thicknessProfile, differencing=differencing, solver=solver, nOutputs=nOutputs, terminateAt=terminateAt)
println("finished sim")

#%%

stoppoint = -20

midpoint = (length(sol1.u)+stoppoint)÷2
C_peak, ind_peak = findmax(reshape(sol1.u[midpoint], p1.dims...)[:,dims[2]÷2])
νs   = collect(range(0.0, 1.0, p1.dims[1]))
ν_peak = νs[ind_peak]
Ẽ = p1.K₂/(1+p1.K₂)
D = Ẽ*p1.K₂*K₄/(1+α_C)
t̃₀ = sol1.t[midpoint] - 1/(4.0*π*D*C_peak^2)
ν₀ = ν_peak - Ẽ*β*(sol1.t[midpoint]-t̃₀)/(1+α_C)
νsOffset = νs.-ν₀
tsOffset = sol1.t.-t̃₀
firstPositivetIndex = findfirst(x->x>0, tsOffset)

midpoint = length(sol1.u)÷2
C_peak, ind_peak = findmax(reshape(sol1.u[midpoint], p1.dims...)[:,1])
νs   = collect(range(0.0, 1.0, p1.dims[1]))
ν_peak = νs[ind_peak]
Ẽ = p1.K₂/(1+p1.K₂)
D = Ẽ*p1.K₂*K₄/(1+α_C)

t̃₀ = sol1.t[midpoint] - 1/(4.0*π*D*C_peak^2)
ν₀ = ν_peak - Ẽ*β*(sol1.t[midpoint]-t̃₀)/(1+α_C)
νsOffset = νs.-ν₀
tsOffset = sol1.t.-t̃₀
firstPositivetIndex = findfirst(x->x>0, tsOffset)

midpoint = (length(sol1.u)+stoppoint)÷2
C_peak, ind_peak = findmax(reshape(sol1.u[midpoint], p1.dims...)[:,1])
νs   = collect(range(0.0, 1.0, p1.dims[1]))
ν_peak = νs[ind_peak]
Ẽ = p1.K₂/(1+p1.K₂)
D = Ẽ*p1.K₂*K₄/(1+α_C)

t̃₀ = sol1.t[midpoint] - 1/(4.0*π*D*C_peak^2)
ν₀ = ν_peak - Ẽ*β*(sol1.t[midpoint]-t̃₀)/(1+α_C)
νsOffset = νs.-ν₀
tsOffset = sol1.t.-t̃₀
firstPositivetIndex = findfirst(x->x>0, tsOffset)

# fig = Figure(size=(1000,1000), fontsize=32)
# ax1 = CairoMakie.Axis(fig[1, 1])
# ax1.xlabel = L"\nu"
# ax1.ylabel = L"\tilde{C}"
# analyticLine = Observable(zeros(dims[1]))
# numericLine = Observable(zeros(dims[1]))
# l1 = lines!(ax1, νs, analyticLine, color=(:red, 0.75), linewidth=4)
# l2 = lines!(ax1, νs, numericLine, color=(:blue, 0.75), linewidth=4, linestyle=:dot)
# Legend(fig[1,2], [l1, l2], ["Asymptotic", "Numeric"])
# ylims!(ax1, (-2.0, maximum(sol1.u[1])))
# xlims!(ax1, (0.0, 1.0))
# # analyticVals = zeros(size(νsOffset)) 

# ax2 = CairoMakie.Axis(fig[2, 1])#, aspect=1)
# ax2.xlabel = L"\tilde{t}"
# ax2.ylabel = L"\tilde{M}_\phi"
# xlims!(ax2, (0.0, sol1.t[end]))
# ylims!(ax2, (0.0, 1.05*π))
# nFrames = 101
# analyticMs = Observable(zeros(nFrames))
# numericMs = Observable(zeros(nFrames))
# # Ms3 = Observable(zeros(nFrames))
# Ts = Observable(zeros(nFrames))
# allLines = [ lines!(ax2, Ts, analyticMs, color=(:red, 0.75), linewidth=4), 
#                 lines!(ax2, Ts, numericMs, color=(:blue, 0.75), linewidth=4, linestyle=:dot) ]
# allLabels = ["Asymptotic", "Numeric"]
# # l5 = lines!(ax2, Ts, Ms3, color=(:green, 0.75), linewidth=4, linestyle=:dot)
# Legend(fig[2,2], allLines, allLabels)
# record(fig, datadir("sims",subFolder, folderName, "analyticCs2.mp4"), 1:length(sol1.t)÷nFrames; framerate=10) do i
#     if tsOffset[i] > 0
#         # analyticVals .= homogeneousWidthC.(νsOffset, K₂, K₄, α_C, β, tsOffset[(i-1)*nFrames+1])
#         analyticLine[] .= homogeneousWidthC.(νsOffset, K₂, K₄, α_C, β, tsOffset[(i-1)*nFrames+1]) # analyticVals
#         uInternal = reshape(sol1.u[(i-1)*nFrames+1], dims...)
#         numericLine[] .= uInternal[:,dims[2]÷2]
#         analyticLine[] = analyticLine[]
#         numericLine[] = numericLine[]

#         Ts[][i] = sol1.t[(i-1)*nFrames+1]
#         Ts[] = Ts[]
#         analyticMs[][i] = M̃ϕAnalytic(ϕ, ν₀, sol1.t[(i-1)*nFrames+1]-t̃₀, α_C, β, p1.K₂, K₄)
#         numericMs[][i] = M̃ϕ(sol1.u[(i-1)*nFrames+1], p1.W, p1.dims, p1.dν, p1.hᵥ, ϕ) 
#         analyticMs[] = analyticMs[]
#         numericMs[] = numericMs[]
#     end
# end


#%%

fig = Figure(size=(1000,800), fontsize=18, figure_padding=25)

g1 = GridLayout(fig[1,1])
g2 = GridLayout(fig[2,1])

ax1 = Axis(g1[1, 1])
# Label(g1[2,1, Top()], "(a)")
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
# axislegend(ax1, allLines_ax1, labels_ax1, labelsize = 16)
ax1.xlabel = L"\nu"
ax1.ylabel = L"\tilde{C}"
# Label(g1[1, 1, Bottom()], L"\nu")
# Label(g1[1, 1, Left()], L"\tilde{C}", rotation=π/2)
ax1.yticks = (0.0:20.0:20.0, [L"0.0", L"20.0"])
ax1.xticks = (0.0:1.0:1.0, [L"0.0", L"1.0"])
text!(ax1, Point{2,Float64}(0.95,0.9*20.0), text="(a)", color=:black, align=(:right, :bottom)) 
text!(ax1, Point{2,Float64}(0.05, 0.9*20.0), text = L"t=%$(allTs_ax1[1])", color=:red) 
text!(ax1, Point{2,Float64}(0.35, 0.4*20.0), text = L"t=%$(allTs_ax1[2])", color=:green) 
text!(ax1, Point{2,Float64}(0.7, 0.3*20.0), text = L"t=%$(allTs_ax1[3])", color=:blue) 
ylims!(ax1, (0.0, 20.0))
xlims!(ax1, (0.0, 1.0))




ax2 = Axis(g1[1,2], yticks = (0.0:π/2.0:π, [L"0", L"π/2", L"π"]))
# Label(g1[2,2, Top()], "(b)")
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
l5 = vlines!(ax2, sol1.t[ind], color=(:black, 0.5))#, linewidth=4
tEndString = @sprintf("%.2f", sol1.t[end+stoppoint])
ax2.xticks = ([0.0, sol1.t[ind], sol1.t[end+stoppoint]], [L"0.0", L"\tilde{T}_{r50}", L"%$(tEndString)"])
ax2.xlabel = L"\tilde{t}"
ax2.ylabel = L"\tilde{M}_\phi"
text!(ax2, Point{2,Float64}(0.95*0.29,0.9*π), text="(b)", color=:black, align=(:right, :bottom)) 

# ax2 = Axis(g1[1,2], yticks = (0.0:π/2.0:π, [L"0", L"π/2", L"π"]))
# Label(g1[2,2, Top()], "(b)")
# ylims!(ax2, (0.0, 1.05*π))
# xlims!(ax2, (0.0, sol1.t[end+stoppoint]))
# tSeries = sol1.t[firstPositivetIndex:end+stoppoint]
# numericalMs1floor = [M̃ϕ(u, p1.W, p1.dims, p1.dν, p1.hᵥ, ϕ, thresh="floor") for u in sol1.u[firstPositivetIndex:end+stoppoint]]
# numericalMs1ceil = [M̃ϕ(u, p1.W, p1.dims, p1.dν, p1.hᵥ, ϕ, thresh="ceil") for u in sol1.u[firstPositivetIndex:end+stoppoint]]
# numericalMs2floor = [M̃ϕ2(u, p1.W, p1.dims, p1.dν, p1.hᵥ, ϕ, thresh="floor") for u in sol1.u[firstPositivetIndex:end+stoppoint]]
# numericalMs2ceil = [M̃ϕ2(u, p1.W, p1.dims, p1.dν, p1.hᵥ, ϕ, thresh="ceil") for u in sol1.u[firstPositivetIndex:end+stoppoint]]
# # analyticMs = [M̃ϕAnalytic.(ϕ, ν₀, τ, α_C, β, p1.K₂, K₄) for τ in tsOffset[firstPositivetIndex:end+stoppoint]]
# allLines = [ lines!(ax2, tSeries, numericalMs1floor, linewidth=4, color=(:red, 0.5)), 
#                 lines!(ax2, tSeries, numericalMs1ceil, linewidth=4, color=(:green, 0.5)), 
#                 lines!(ax2, tSeries, numericalMs2floor, linewidth=4, color=(:blue, 0.5)), 
#                 lines!(ax2, tSeries, numericalMs2ceil, linewidth=4, color=(:black, 0.5)),
#             ]
# allLabels = [ "numericalMs1floor",
#     "numericalMs1ceil",
#     "numericalMs2floor",
#     "numericalMs2ceil",
# ]
# # l4 = lines!(ax2, tSeries, analyticMs, linestyle=:dot , linewidth=4, color=(:red, 1.0))
# ind = findfirst(x->M̃ϕ(x, p1.W, p1.dims, p1.dν, p1.hᵥ, ϕ)>=π/2.0, sol1.u)
# l5 = vlines!(ax2, sol1.t[ind], color=(:black, 0.5))#, linewidth=4
# tEndString = @sprintf("%.2f", sol1.t[end+stoppoint])
# ax2.xticks = ([0.0, sol1.t[ind], sol1.t[end+stoppoint]], [L"0.0", L"\tilde{T}_{r50}", L"%$(tEndString)"])
# ax2.xlabel = L"\tilde{t}"
# ax2.ylabel = L"\tilde{M}_\phi"
# axislegend(ax2, allLines, allLabels, labelsize = 16, position = :lt)


linesVec_ax3 = []
labelsVec_ax3 = []
ax3 = Axis(g2[1,1])
linesVec_ax3 = []
labelsVec_ax3 = []
hcutoff = (2.0*k_Sa/k_Sd)*((𝒮*k₁*k₃)/(2.0*Ωperp*k₂*k₄) - 1.0)
push!(linesVec_ax3, lines!(ax3, h₀s, 𝒫sim, color=(:red, 0.5), linewidth=4))
push!(labelsVec_ax3, "Numeric")
push!(linesVec_ax3, lines!(ax3, [0.0, h₀s..., hcutoff], 𝒫analytic, color=(:blue, 1.0), linewidth=4, linestyle=:dot))
push!(labelsVec_ax3, "Asymptotic")
push!(linesVec_ax3, vlines!(ax3, h_C, color=(:black, 0.5)))#, linewidth=4))
push!(linesVec_ax3, vlines!(ax3, h_S, color=(:black, 0.5)))#, linewidth=4))
push!(linesVec_ax3, vlines!(ax3, hcutoff, color=(:black, 0.5)))#, linewidth=4))

ax3.xticks = ([0.0, h_C, h_S, hcutoff], [L"0", L"  h_C", L"h_S", L"h_{cut-off}"])
ax3.yticks = ([0.0, 0.0001, 0.0002, 0.0003, 0.0004], [L"0.0", L"1.0", L"2.0", L"3.0", L"4.0"])
ax3.xaxis.elements[:ticklabels].align = tuple.([:right, :left, :center, :center], :top)
xlims!(ax3, (0.0, 1.05*hcutoff))
ylims!(ax3, (0.0, 1.1*maximum(𝒫analytic)))
ax3.xlabel = L"h_0"
ax3.ylabel = L"𝓟^*_{50}/10^{-4}"

text!(ax3, Point{2,Float64}(0.9*1.05*hcutoff, 0.9*1.1*maximum(𝒫analytic)), text="(c)", color=:black, align=(:right, :bottom)) 

colsize!(g2, 1, Aspect(1, 1.5))
resize_to_layout!(fig)

save(datadir("sims", subFolder, folderName, "Figure2_ν0=$(@sprintf("%.6f", ν₀))_t0=$(@sprintf("%.6f", t̃₀)).png"), fig)
display(fig)

@show t̃₀
@show ν₀
