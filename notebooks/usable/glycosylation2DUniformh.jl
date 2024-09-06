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

# δC = (1 + α_C*hᵥ)⁻¹ * ( E*K_2*K_4*C_νν  -  E*β*C_ν  + D*∇_perp⋅(h*∇_perp*C) )
# - E*β*C_ν  --> E*∇cdot*β*Pν*Aᵤₚ*Cᵥ                β*Aᵤₚ*Cᵥ gives the velocity of each edge multiplied by the upstream concentration of each edge, giving flux on edge. Adding Pν excludes x-directed edges.
# + E*K_2*K_4*C_νν --> E*K_2*K_4*∇cdot*Pν*∇ₑ*Cᵥ      This is a diffusion term with diffusion constant K_2*K_4, which does not vary over space. ∇ₑ*C gives diffusive flux on each edge; Pν picks only ν directed edges, ∇cdot takes divergence over those edges
# + D * ∇_perp ⋅ (h ∇_perp * C ) --> ∇cdot*hₑ*Px*𝓓ₑ*∇ₑ*Cᵥ  Another diffusive term. 𝓓ₑ*∇ₑ*Cᵥ gives diffusive flux over all edges; Px picks out only x-directed edges; hₑ multiplies by thickness at each edge; ∇cdot takes divergence over those edges


#


using OrdinaryDiffEq
using SparseArrays
using UnPack
using CairoMakie 
using FromFile
using DrWatson
using Printf
using SciMLOperators
using DataFrames
using Dates
using LinearAlgebra

@from "$(srcdir("MakeIncidenceMatrix.jl"))" using MakeIncidenceMatrix
@from "$(srcdir("MakeWeightMatrices.jl"))" using MakeWeightMatrices
@from "$(srcdir("Visualise.jl"))" using Visualise
@from "$(srcdir("UsefulFunctions.jl"))" using UsefulFunctions

Ωperp = 100.0  # Lumen footprint area
N     = 100         # Maximum polymer length 
k_Cd = 0.9 # Complex desorption rate
k_Ca = 1.1 # Complex adsorption rate
k_Sd = 0.9 # Substrate desorption rate
k_Sa = 1.1 # Substrate adsorption rate
k₁   = 1.0   # Complex formation forward reaction rate 
k₂   = 0.6   # Complex dissociation reverse reaction rate 
k₃   = 1.1   # Product formation
k₄   = 0.6  # Product dissociation 
E_0 = 1.0
𝓒 = 100.0
𝓢 = 1000.0
D_C  = 1.0  # Monomer/polymer diffusivity
D_S  = 1.0  # Substrate diffusivity
Tᵣstar  = 100.0  # Release time

λ = (𝓢/(2*Ωperp))*(k₁*k₃/(k₂*k₄))
@show λ

𝓔    = 2*Ωperp*E_0   # Total enzyme mass
K₃  = k₃/k₁    # Non-dimensionalised product formation rate
K₄  = k₄/k₁    # Non-dimensionalised prodict dissociation rate
δ_C = π*D_C/(k₁*𝓔)
δ_S = π*D_S/(k₁*𝓔)
Tᵣ  = k₁*𝓔*Tᵣstar/(2*Ωperp)
ϕ = 0.5

# PDE discretisation parameters 
Nx       = 101             # Number of discretisation points in space
Nν       = 101             # Number of discretisation points in polymerisation space
Nghost   = 1           # Number of ghost points on each side of the domain 
Nνplus   = Nν+2*Nghost # Number of discretised points including ghost points 
Nxplus   = Nx+2*Nghost # Number of discretised points including ghost points 
dimsPlus = [Nνplus, Nxplus]
dims     = (Nν, Nx)
# xMax     = 100.0
xs       = collect(range(0.0, Ωperp, dimsPlus[2])) # Positions of discretised vertices in space
dx       = xs[2]-xs[1]
νMax     = 1.0
νs       = collect(range(0.0, νMax, Nνplus)) # Positions of discretised vertices in polymerisation space 
dν       = νs[2]-νs[1]
spacing  = (dν, dx)

# Incidence matrices 
A   = makeIncidenceMatrix3D(Nνplus, Nxplus, 1)
Ā   = abs.(A)
Aᵀ  = transpose(A)
Āᵀ  = transpose(Ā)
# Change of concentration at a given vertex due to advection is the divergence of advective flux on each adjacent edge.
# When calculting advective flux into and out of a vertex k, we need a value of concentration for each edge j entering or leaving that vertex.
# We could use the mean concentration of its adjacent vertices 0.5.*Ā*Cᵥ, but it is better to take the concentration of the upstream vertex.
# This means mapping the value of the vertex k for which A[j,k]=-1, meaning the vertex that edge j exits, to edge j.
# To pick out these values we use Aᵤₚ*Cᵥ where Aᵤₚ=(Ā-A)/2
Aᵤₚ = dropzeros((Ā-A).÷2)   

# Number of vertices and number of edges, total and in each dimension
nVerts  = Nνplus*Nxplus       # Total number of vertices 
nEdgesi = (Nνplus-1)*Nxplus  # Number of i-directed edges (ν, in this case)
nEdgesj = Nνplus*(Nxplus-1)  # Number of j-directed edges (x, in this case)
nEdges  = nEdgesi+nEdgesj     # Total number of edges over all dimensions 

# Ghost point masks; vectors and sparse diagonal matrices to exclude ghost points and edges connected to ghost points 
ghostVertexMaskVec    = makeGhostVertexMask(dimsPlus)
ghostVertexMaskSparse = spdiagm(ghostVertexMaskVec)
ghostEdgeMaskVec      = makeGhostEdgeMask(dimsPlus)
ghostEdgeMaskSparse   = spdiagm(ghostEdgeMaskVec)

# Matrices for picking out ν and xy directions in derivatives 
Pν  = ghostEdgeMaskSparse*spdiagm(vcat(ones(Int64, nEdgesi), zeros(Int64, nEdgesj)))   # Diagonal sparse matrix to exclude all xy edges #########and ν edges adjacent to ghost points  
Px  = ghostEdgeMaskSparse*spdiagm(vcat(zeros(Int64, nEdgesi), ones(Int64, nEdgesj)))   # Diagonal sparse matrix to exclude all ν edges #########and xy edges adjacent to ghost points 

# Weights
W   = vertexVolumeWeightsMatrix(dimsPlus, spacing)
W⁻¹ = vertexVolumeWeightsInverseMatrix(dimsPlus, spacing)
l⁻¹ = edgeLengthInverseMatrix(dimsPlus, spacing) # Diagonal matrix of edge lengths

# Gradient operators 
∇ₑ = l⁻¹*A       # Gradient operator giving gradient on each edge
∇cdot = -W⁻¹*Aᵀ  # Divergence operator giving divergence on each vertex calculated from edges 

# Diagonal matrix of areas perpendicular to each edge, 
# meaning the area through which diffusive flux in the direction of a given edge passes
Aperpₑ = edgePerpendicularAreaMatrix(dimsPlus, spacing)

# Set value of Fₑ at each point in space
matFₑ = ones(Float64, Nxplus)
integF = dx*sum(matFₑ[2:end-1])
# Ensure integral of Fₑ over space is π
matFₑ .*= π/integF
matE = zeros(Nxplus)
Esparse = spzeros(nVerts, nVerts)

h₀s = collect(0.1:0.1:3.0)

sols = []
hᵥs = []
α_Cs = []
Ωs =[]
C_bs =[]

for h₀ in h₀s
    @show h₀

    Ω     = h₀*Ωperp         # Lumen volume
    α_C = (k_Cd*Ω)/(2*k_Ca*Ωperp) # Balance of complex in bulk to complex on membrane       units of m²?
    α_S = (k_Sd*Ω)/(2*k_Sa*Ωperp) # Balance of substrate in bulk to substrate on membrane   units of m²?
    C_b  = 𝓒/Ω 
    S_b  = 𝓢/Ω 
    C_0 = C_b*h₀/(2*(1+α_C))      # Early surface monomer concentration
    S_0 = S_b*h₀/(2*(1+α_S))      # Early surface substrate concentration 
    K₂  = k₂/(k₁*C_0)              # (k₂/(k₁*C_b))*((2*k_Ca*Ωperp + k_Cd*Ω)/(k_Ca*Ω)) # Non-dimensionalised complex formation net reaction rate
    σ   = S_0/C_0                         #(k_Sa*S_b*(2*k_Ca*Ωperp + k_Cd*Ω)) / (k_Ca*C_b*(2*k_Sa*Ωperp + k_Sd*Ω))
    ϵ   = E_0/C_0                  # 𝓔*(2*k_Ca*Ωperp + k_Cd*Ω) / (2*k_Ca*C_b*Ωperp)
    𝓓   = α_C*δ_C*N^2*(K₂ + σ*K₃)
    β = N*(σ*K₃ - K₂*K₄)

    𝓓ₑ       = 𝓓.*ghostEdgeMaskSparse*Aperpₑ # Sparse diagonal matrix of diffusivities over edges 
    
    mat_h    = h₀.*ones(dimsPlus...)
    # Diagonal matrices of compartment thickness h over all vertices hᵥ
    # Also diagonal matrix of thickness over edges, formed by taking mean of h at adjacent vertices 0.5.*Ā*hᵥ
    hᵥ_vec = reshape(mat_h, nVerts)         # Cisternal thickness evaluated over vertices 
    hₑ_vec = 0.5.*Ā*hᵥ_vec                  # Cisternal thickness evaluated over edges (mean of adjacent vertices)
    hᵥ = spdiagm(hᵥ_vec)                    # Cisternal thickness over vertices, as a sparse diagonal matrix
    hₑ = spdiagm(hₑ_vec)                    # Cisternal thickness over edges, as a sparse diagonal matrix
    aᵥ = spdiagm(1.0./(1.0 .+ α_C.*hᵥ_vec)) # Prefactor 1/(1+α_C*hᵥ(x)) evaluated over vertices, packaged into a sparse diagonal matrix for convenience
    aₑ = spdiagm(1.0./(1.0 .+ α_C.*hₑ_vec)) # Prefactor 1/(1+α_C*hₑ(x)) evaluated over edges, packaged into a sparse diagonal matrix for convenience

    push!(hᵥs, hᵥ)
    push!(Ωs, h₀*Ωperp)
    push!(α_Cs, α_C)
    push!(C_bs, C_b)

    # Initial conditions using Gaussian
    # Assume initially that we only have cargo in bulk, so normalisation with 𝓒 is done only using bulk concentration
    uMat = zeros(Float64, Nνplus, Nxplus)
    for xx=1:Nxplus, νν=1:Nνplus
        uMat[νν, xx] =  exp(-(νs[νν]^2)/(0.01*νMax)^2)
    end
    u0 = reshape(uMat, nVerts)
    u0[ghostVertexMaskVec.!=true] .= 0.0
    # For integration to normalise the number of monomers, we need to multiply the concentration at each point by the ν value of that point
    νMat = ones(dimsPlus...)
    for ii=1:Nνplus
        νMat[ii,:].*=(ii-1)
    end
    νSparse = spdiagm(reshape(νMat, nVerts))
    integ = sum(ghostVertexMaskSparse*W*νSparse*u0)
    u0 .*= 𝓒/integ

    E!(u0, dimsPlus, Esparse, matE, matFₑ, K₂, dν)

    # PDE operator components
    L1 = aᵥ*∇cdot*(K₂*K₄.*Pν*∇ₑ - β.*Pν*Aᵤₚ)
    L2 = aᵥ*∇cdot*(hₑ*Px*𝓓ₑ*∇ₑ)
    p = (L1=L1, L2=L2, u0=u0, dimsPlus=dimsPlus, Esparse=Esparse, matE=matE, matFₑ=matFₑ, K₂=K₂, dν=dν)
    L = MatrixOperator(Esparse*L1.+L2, update_func! = updateOperator!)
    prob = ODEProblem(L, u0, (0.0, Tᵣ), p)
    sol = solve(prob, Vern9(), saveat=(Tᵣ)/100.0)

    push!(sols, sol.u)
end

#%%

fig = Figure(size=(500,500))
ax = Axis(fig[1,1])
Pstars = Float64[]
for i=1:length(sols)
    push!(Pstars, P_star(sols[i][end], W, ghostVertexMaskVec, dims, hᵥs[i], ϕ, α_Cs[i], C_bs[i], Ωs[i], Tᵣstar))
end
lines!(ax, h₀s, Pstars)
ax.xlabel = "h₀"
ax.ylabel = L"𝓟^*"
display(fig)
save("simulationPvsh.png",fig)

# fig = Figure(size=(500,500))
# ax = Axis(fig[1,1])
# Mstars = Float64[]
# for i=1:length(uFinals)
#     push!(Mstars, M_star(uFinals[i], W, ghostVertexMaskVec, dims, hᵥs[i], ϕ, α_Cs[i], C_b, Ωs[i]))
# end
# lines!(ax, h₀s, Mstars)
# ax.xlabel = "h₀"
# ax.ylabel = L"M^*"
# display(fig)
# save("simulationMvsh.png",fig)

# fig = Figure(size=(1000,1000))
# ax  = Axis3(fig[1, 1], aspect=:equal, azimuth=-π/4)
# ax.xlabel = "ν"
# ax.ylabel = "x"
# ax.zlabel = "c"
# for (i,u) in enumerate(uFinals)
#     empty!(ax)
#     uInternal = reshape(uFinals[1][ghostVertexMaskVec], dims)
#     surface!(ax, xs[Nghost+1:end-Nghost], νs[Nghost+1:end-Nghost], uInternal, colormap=:batlow)
#     save("uFinalSurface$i.png", fig)
# end
# display(fig)


for j=1:length(sols)
    fig = Figure(size=(1000,1000))
    ax = Axis3(fig[1, 1], aspect=:equal, azimuth=-π/4)
    ax.xlabel = "ν"
    ax.ylabel = "x"
    ax.zlabel = "c"
    uInternal = Observable(zeros(dims))
    globalmin = minimum([minimum(u[ghostVertexMaskVec]) for u in sols[j]])
    globalmax = maximum([maximum(u[ghostVertexMaskVec]) for u in sols[j]])
    zlims!(ax, (globalmin, globalmax))
    clims = (globalmin,globalmax)
    surface!(ax, xs[Nghost+1:end-Nghost], νs[Nghost+1:end-Nghost], uInternal, colorrange=clims, colormap=:batlow)
    record(fig, datadir("concentrationSurfaceMovie$j.mp4"), 1:length(sols[j]); framerate=10) do i
        uInternal[] .= reshape(sols[j][i][ghostVertexMaskVec], dims)
        uInternal[] = uInternal[]
    end
end
# concentrationSurfaceMovie(sol.u, sol.t, xs, νs, (Nν,Nx), Nghost, ghostVertexMaskVec; subFolder=subFolder, folderName=folderName)

# function spaceIntegralOver_ν_Movie(solu, ts, xs, νs, dims, Nghost, vertexWeightsMatrix, ghostVertexMaskVec; subFolder="", folderName="")
#     isdir(datadir("sims", subFolder, folderName)) ? nothing : mkdir(datadir("sims", subFolder, folderName))
#     # Find limits
#     uInternal2D = reshape((vertexWeightsMatrix*solu[end])[ghostVertexMaskVec], dims)
#     M = sum(uInternal2D, dims=2)[:,1]
#     minima = Float64[]
#     maxima = Float64[]
#     for i=1:length(ts)
#         uInternal2D .= reshape((vertexWeightsMatrix*solu[i])[ghostVertexMaskVec], dims)
#         M .= sum(uInternal2D, dims=2)[:,1]
#         push!(minima, minimum(M))
#         push!(maxima, maximum(M))
#     end
#     globalmin = minimum(minima)
#     globalmax = maximum(maxima)

#     fig = Figure(size=(1000,1000))
#     ax = CairoMakie.Axis(fig[1, 1], aspect=1)
#     ax.xlabel = "ν"
#     ax.ylabel = "M, ∱cdxdy"
#     ax.title = "Integral of Cₛ over x against ν"
#     M = Observable(zeros(dims[1]))
#     lines!(ax, νs[1:Nghost:end-2*Nghost], M)
#     ylims!(ax, (globalmin, globalmax))
#     record(fig, datadir("sims",subFolder, folderName, "spaceIntegralOver_ν_Movie.mp4"), 1:length(ts); framerate=10) do i
#         uInternal2D .= reshape((vertexWeightsMatrix*solu[i])[ghostVertexMaskVec], dims)
#         M[] .= sum(uInternal2D, dims=2)[:,1]
#         M[] = M[]
#     end
#     save(datadir("sims",subFolder,folderName,"finalSpaceIntegralOver_ν.png"), fig)
#     return nothing
# end



# # K₁ = 1.0
# # K₂ = 1.0
# # K₃ = 2.0
# # K₄ = 1.0  
# # α_C = 10.0
# # δ_C = 1.0
# # σ = 1.0
# # N = 100
# # β = N*(σ*K₃ - K₂*K₄)
# # 𝓓 = α_C*δ_C*N^2*(K₂+σ*K₃)
# # Tᵣ = 1.0


# # @show K₂ 
# # @show K₃ 
# # @show K₄ 
# # @show α_C 
# # @show δ_C 
# # @show σ 
# # @show N 
# # @show β 
# # @show 𝓓 
# # @show Tᵣ 