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

textheme = Theme(fonts=(; regular=texfont(:text),
                        bold=texfont(:bold),
                        italic=texfont(:italic),
                        bold_italic=texfont(:bolditalic)))

@from "$(srcdir("Glycosylation.jl"))" using Glycosylation
@from "$(srcdir("Visualise.jl"))" using Visualise
@from "$(srcdir("UsefulFunctions.jl"))" using UsefulFunctions
@from "$(srcdir("MakeWeightMatrices.jl"))" using MakeWeightMatrices
@from "$(srcdir("DerivedParameters.jl"))" using DerivedParameters

#%%

subFolder = "Figure2"
# Create directory for run data labelled with current time.
folderName = "$(Dates.format(Dates.now(),"yy-mm-dd-HH-MM-SS"))"
# Create frames subdirectory to store system state at each output time
mkpath(datadir("sims",subFolder,folderName))
terminateAt = "nuWall"
thicknessProfile = "uniform"
differencing = "centre"
solver = SSPRK432()
nOutputs = 10000
σGRF = 0.2

nSpatialDims = 1
Ngrid = 201
dims = fill(Ngrid, nSpatialDims+1)

include(projectdir("notebooks", "paramsRaw.jl"))

#%%

h_C = 2*k_Ca/k_Cd
h_S = 2*k_Sa/k_Sd
hMax = h_C*5
hMin = h_C/10
h₀s = collect(hMin:2*hMin:hMax)

h₀s2 = collect(hMax+1.0:1.0:7.0)
append!(h₀s, h₀s2)

Ωs = h₀s.*Ωperp      # Dimensional lumen volume 

#%%

𝒫sim = []
𝒫simEq50 = []
𝒫analytic = []
sols = []
ps = []
for i=1:length(h₀s)
    @show h₀s[i]    
    derivedParams = derivedParameters(Ωs[i], Ωperp, N, k_Cd, k_Ca, k_Sd, k_Sa, k₁, k₂, k₃, k₄, 𝒞, 𝒮, ℰ, D_C, D_S, Tᵣstar; checks=false)
    @unpack L₀, E₀, C_b, S_b, δ_C, δ_S, α_C, α_S, C₀, S₀, Tᵣ, T̃ᵣ, K₂, K₃, K₄, σ, ϵ, 𝒟, β, h_C, h_S, λ, ζ, γ, Δ, F = derivedParams
    sol, p = glycosylationAnyD(dims, K₂, K₄, 1000.0, α_C, 𝒟, β, thickness=thicknessProfile, differencing=differencing, solver=solver, nOutputs=nOutputs, terminateAt="halfProduction", saveIntermediate=false) 
    Tᵣ₅₀Star = sol.t[end]*(N^2)*(K₂+σ*K₃)/(k₁*E₀)
    push!(sols, sol)
    push!(ps, p)
    push!(𝒫sim, Mstarϕ(sol.u[end], p.W, p.dims, p.dν, p.hᵥ, α_C, 𝒞, ϕ)/Tᵣ₅₀Star)
    push!(𝒫analytic, 𝒫star₅₀Analytic(h₀s[i], h_C, h_S, k₁, k₂, k₃, k₄, Ωperp, 𝒮, 𝒞, ℰ, N, ϕ))
    push!(𝒫simEq50, 𝒫star₅₀Numeric(N, k₁, k₂, k₃, 𝒞, ℰ, 𝒮, h₀s[i], k_Ca, k_Cd, k_Sa, k_Sd, Ωperp, sol.t[end]))    
end

#%%

derivedParams = derivedParameters(Ω, Ωperp, N, k_Cd, k_Ca, k_Sd, k_Sa, k₁, k₂, k₃, k₄, 𝒞, 𝒮, ℰ, D_C, D_S, Tᵣstar; checks=false)
@unpack L₀, E₀, C_b, S_b, δ_C, δ_S, α_C, α_S, C₀, S₀, Tᵣ, T̃ᵣ, K₂, K₃, K₄, σ, ϵ, 𝒟, β, h_C, h_S, λ, ζ, γ, Δ, F = derivedParams

# #%%

sol1, p1 = glycosylationAnyD(dims, K₂, K₄, T̃ᵣ, α_C, 𝒟, β, thickness=thicknessProfile, differencing=differencing, solver=solver, nOutputs=nOutputs, terminateAt=terminateAt)
println("finished sim")


#%%


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

fig = Figure(size=(1000,1000), fontsize=32, theme=textheme)
ax1 = CairoMakie.Axis(fig[1, 1])
ax1.xlabel = L"\nu"
ax1.ylabel = L"\tilde{C}"
analyticLine = Observable(zeros(dims[1]))
numericLine = Observable(zeros(dims[1]))
l1 = lines!(ax1, νs, analyticLine, color=(:red, 0.75), linewidth=4)
l2 = lines!(ax1, νs, numericLine, color=(:blue, 0.75), linewidth=4, linestyle=:dot)
Legend(fig[1,2], [l1, l2], ["Analytic", "Numeric"])
ylims!(ax1, (-2.0, maximum(sol1.u[1])))
xlims!(ax1, (0.0, 1.0))
analyticVals = zeros(size(νsOffset)) 

ax2 = CairoMakie.Axis(fig[2, 1])#, aspect=1)
ax2.xlabel = L"\tilde{t}"
ax2.ylabel = L"\tilde{M}_\phi"
xlims!(ax2, (0.0, sol1.t[end]))
ylims!(ax2, (0.0, 1.05*π))
Ms = Observable(zeros(length(sol1.t)))
Ms2 = Observable(zeros(length(sol1.t)))
Ts = Observable(zeros(length(sol1.t)))
l3 = lines!(ax2, Ts, Ms, color=(:red, 0.75), linewidth=4)
l3 = lines!(ax2, Ts, Ms2, color=(:blue, 0.75), linewidth=4, linestyle=:dot)
record(fig, datadir("sims",subFolder, folderName, "analyticCs.mp4"), 1:length(sol1.t)÷101:length(sol1.t); framerate=20) do i
    if tsOffset[i] > 0
        analyticVals .= homogeneousWidthC.(νsOffset, K₂, K₄, α_C, β, tsOffset[i])
        analyticLine[] .= analyticVals
        uInternal = reshape(sol1.u[i], dims...)
        numericLine[] .= uInternal[:,dims[2]÷2]
        analyticLine[] = analyticLine[]
        numericLine[] = numericLine[]

        Ts[][i] = sol1.t[i]
        Ts[] = Ts[]
        Ms[][i] = M̃ϕ(sol1.u[i], p1.W, p1.dims, p1.dν, p1.hᵥ, ϕ)
        Ms2[][i] = M̃ϕAnalytic.(ϕ, ν₀, sol1.t[i]-t̃₀, α_C, β, p1.K₂, K₄)
        Ms[] = Ms[]
        Ms2[] = Ms2[]
    end
end



#%%

fig = Figure(size=(1000,1500), fontsize=32, figure_padding = 30)#, theme=textheme)

ax1 = Axis(fig[1, 1])
Label(fig[2,1, Top()], "(a)")
allLines_ax1 = []
allTs_ax1 = []
colorsUsed = [(:red), (:green), (:blue)]
for (c,i) in enumerate([firstPositivetIndex, (length(sol1.t)-firstPositivetIndex)÷2+firstPositivetIndex, length(sol1.t)])
    uInternal = reshape(sol1.u[i], p1.dims...)
    push!(allLines_ax1, lines!(ax1, νs, uInternal[:,1], linestyle=:solid, color=(colorsUsed[c], 0.5), linewidth=4))
    push!(allLines_ax1, lines!(ax1, νs, homogeneousWidthC.(νsOffset, p1.K₂, K₄, α_C, β, tsOffset[i]), linestyle=:dot, color=(colorsUsed[c], 1.0), linewidth=4))
    push!(allTs_ax1, @sprintf("%.2f", tsOffset[i]))
end
labels_ax1 = []
for t in allTs_ax1
    push!(labels_ax1, "Numeric, t=$t")
    push!(labels_ax1, "Analytic, t=$t")
end
axislegend(ax1, allLines_ax1, labels_ax1, labelsize = 16)
ax1.xlabel = L"\nu"
ax1.ylabel = L"\tilde{C}"
ylims!(ax1, (0.0, 20.0))
xlims!(ax1, (0.0, 1.0))

ax2 = Axis(fig[3,1], yticks = (0.0:π/2.0:π, [L"0", L"π/2", L"π"]))
Label(fig[4,1, Top()], "(b)")
ylims!(ax2, (0.0, 1.05*π))
xlims!(ax2, (0.0, sol1.t[end]))
tSeries = sol1.t[firstPositivetIndex:end]
numericalMs = [M̃ϕ(u, p1.W, p1.dims, p1.dν, p1.hᵥ, ϕ) for u in sol1.u[firstPositivetIndex:end]]
analyticMs = [M̃ϕAnalytic.(ϕ, ν₀, τ, α_C, β, p1.K₂, K₄) for τ in tsOffset[firstPositivetIndex:end]]
l3 = lines!(ax2, tSeries, numericalMs, linewidth=4, color=(:black, 0.5))
l4 = lines!(ax2, tSeries, analyticMs, linestyle=:dot , linewidth=4, color=(:black, 1.0))
ind = findfirst(x->M̃ϕ(x, p1.W, p1.dims, p1.dν, p1.hᵥ, ϕ)>=π/2.0, sol1.u)
l5 = vlines!(ax2, sol1.t[ind], linewidth=4, color=(:black, 0.5))
tEndString = @sprintf("%.2f", sol1.t[end])
ax2.xticks = ([0.0, sol1.t[ind], sol1.t[end]], [L"0.0", L"\tilde{T}_{r50}", L"%$(tEndString)"])
ax2.xlabel = L"\tilde{t}"
ax2.ylabel = L"\tilde{M}_\phi"
axislegend(ax2, [l3, l4], ["Numeric", "Analytic"], labelsize = 16)


linesVec_ax3 = []
labelsVec_ax3 = []

ax3 = Axis(fig[5,1])
Label(fig[6,1, Top()], "(c)")
linesVec_ax3 = []
labelsVec_ax3 = []
push!(linesVec_ax3, lines!(ax3, h₀s, 𝒫sim, color=:blue, linewidth=4))
push!(labelsVec_ax3, "Numeric")
# push!(linesVec_ax3, lines!(ax3, h₀s, 𝒫simEq50, color=:red), linewidth=4)
# push!(labelsVec_ax3, "Equation 50")
push!(linesVec_ax3, lines!(ax3, h₀s, 𝒫analytic, color=:red, linewidth=4))
push!(labelsVec_ax3, "Asymptotic")
push!(linesVec_ax3, vlines!(ax3, h_C, color=(:black, 0.5), linewidth=4))
# push!(labelsVec_ax3, L"h_C")
push!(linesVec_ax3, vlines!(ax3, h_S, color=(:black, 0.5), linewidth=4))
# push!(labelsVec_ax3, L"h_S")
hcutoff = (2.0*k_Sa/k_Sd)*((𝒮*k₁*k₃)/(2.0*Ωperp*k₂*k₄) - 1.0)
push!(linesVec_ax3, vlines!(ax3, hcutoff, color=(:black, 0.5), linewidth=4))
# push!(labelsVec_ax3, L"h_{cut-off}")

ax3.xticks = ([0.0, h_C, h_S, hcutoff], [L"0", L"  h_C", L"h_S", L"h_{cut-off}"])
ax3.yticks = ([0.0, 0.0001, 0.0002, 0.0003, 0.0004], [L"0.0", L"1.0", L"2.0", L"3.0", L"4.0"])

xlims!(ax3, (0.0, 1.1*maximum([maximum(h₀s), h_C, h_S, hcutoff])))
ylims!(ax3, (0.0, 1.1*maximum(𝒫analytic)))

ax3.xlabel = "h₀"
ax3.ylabel = L"𝓟^*_{50}/10^{-4}"

axislegend(ax3, linesVec_ax3[1:2], labelsVec_ax3, labelsize = 16)

colsize!(fig.layout, 1, Aspect(1, 1.5))
rowsize!(fig.layout, 2, Relative(0.01))
rowsize!(fig.layout, 4, Relative(0.01))
rowsize!(fig.layout, 6, Relative(0.01))

resize_to_layout!(fig)
save(datadir("sims", subFolder, folderName, "Figure2_ν0=$(@sprintf("%.6f", ν₀))_t0=$(@sprintf("%.6f", t̃₀)).png"), fig)
display(fig)

@show t̃₀
@show ν₀

#%%

# 𝒫simEq50 = 𝒫star₅₀Numeric.(N, k₁, k₂, k₃, 𝒞, ℰ, 𝒮, h₀s, k_Ca, k_Cd, k_Sa, k_Sd, Ωperp, [s.t[end] for s in sols])

# fig2 = Figure(size=(1000,500), fontsize=32, figure_padding = 30)#, theme=textheme)
# ax21 = Axis(fig2[1, 1])
# # ax22 = Axis(fig2[2, 1])
# # ax23 = Axis(fig2[3, 1])
# l1 = lines!(ax21, h₀s, 𝒫sim, color=(:black, 0.5), linewidth=8)
# # Label(fig2[1, 1, Top()], L"M^*_\phi/T^*_{r50}")
# # l2 = lines!(ax21, h₀s, 𝒫analytic, color=:green, linewidth=8)
# # Label(fig2[2, 1, Top()], L"Equation 57")
# l3 = lines!(ax21, h₀s, 𝒫simEq50, color=(:red, 1.0), linestyle=:dot, linewidth=8)
# # Label(fig2[3, 1, Top()], L"Equation 50")

# axislegend(ax21, [l1, l2, l3], ["Eq46/48", "Eq57", "Eq50"], labelsize = 16)

# display(fig2)
# save(datadir("sims", subFolder, folderName, "Figure2b.png"), fig2)