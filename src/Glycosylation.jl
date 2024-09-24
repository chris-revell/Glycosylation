#
#  Glycosylation.jl
#  Glycosylation
#
#  Created by Christopher Revell on 09/09/2024.


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



module Glycosylation

using OrdinaryDiffEq
using SparseArrays
using UnPack
using FromFile
using DrWatson
using SciMLOperators
using DataFrames
using Statistics
using InvertedIndices

@from "$(srcdir("MakeIncidenceMatrix.jl"))" using MakeIncidenceMatrix
@from "$(srcdir("MakeWeightMatrices.jl"))" using MakeWeightMatrices
# @from "$(srcdir("Visualise.jl"))" using Visualise
@from "$(srcdir("UsefulFunctions.jl"))" using UsefulFunctions
# @from "$(srcdir("CisternaWidth.jl"))" using CisternaWidth
@from "$(srcdir("DerivedParameterChecks.jl"))" using DerivedParameterChecks


function glycosylationAnyD(xs, mat_h, nSpatialDims, Ngrid, Nghost, Ωperp, N, k_Cd, k_Ca, k_Sd, k_Sa, k₁, k₂, k₃, k₄, E_0, 𝓒, 𝓢, D_C, D_S, Tᵣstar, ϕ)

    
# PDE discretisation parameters 
    Nνplus   = Ngrid
    Nxplus   = Ngrid
    Nyplus   = Ngrid

    nSpatialDims == 1 ? dimsPlus = [Nνplus, Nxplus] : dimsPlus = [Nνplus, Nxplus, Nyplus]
    nSpatialDims == 1 ? dimsReal = [Ngrid, Ngrid] : dimsReal = [Ngrid, Ngrid, Ngrid]

    # Generate xMax and width profile from data or function 
    # xMax, mat_h = hFromData(dimsPlus; cisternaSeriesID=1)
    xMax = maximum(xs)
    dx   = xs[2]-xs[1]
    if nSpatialDims > 1 
        yMax = xMax
        ys   = collect(range(0.0, yMax, Nyplus))
        dy   = ys[2]-ys[1]
    end
    νMax = 1.0
    νs   = collect(range(0.0, νMax, Nνplus)) # Positions of discretised vertices in polymerisation space 
    dν   = νs[2]-νs[1]

    nSpatialDims == 1 ? spacing  = [dν, dx] : spacing  = [dν, dx, dy]

    h₀ = mean(selectdim(mat_h, 1, 1))

    derivedParams = derivedParameterNoChecks(h₀, Ωperp, N, k_Cd, k_Ca, k_Sd, k_Sa, k₁, k₂, k₃, k₄, E_0, 𝓒, 𝓢, D_C, D_S, Tᵣstar; checks=false)
    @unpack 𝓔, K₃, K₄, δ_C, δ_S, Tᵣ, Ω, α_C, α_S, C_b, S_b, C_0, S_0, K₂, σ, ϵ, 𝓓, β, K₂, L₀ = derivedParams 

    λ = (𝓢/(2*Ωperp))*(k₁*k₃/(k₂*k₄))
    @show λ

    A   = makeIncidenceMatrix3D(dimsPlus)
    Ā   = abs.(A)
    Aᵀ  = transpose(A)
    Aᵤₚ = dropzeros((Ā-A).÷2)

    nVerts  = prod(dimsPlus)      # Total number of vertices 
    dimEdgeCount = Int64[]
    for i=1:length(dimsPlus)
        push!(dimEdgeCount, (dimsPlus[i]-1)*prod(dimsPlus[Not(i)]))
    end
    nEdges  = sum(dimEdgeCount)   # Total number of edges over all dimensions 

    # Ghost point masks
    # ghostVertexMaskVec = makeGhostVertexMask(dimsPlus)
    # ghostVertexMaskSparse = spdiagm(ghostVertexMaskVec)
    # ghostEdgeMaskVec = makeGhostEdgeMask(dimsPlus)
    # ghostEdgeMaskSparse = spdiagm(ghostEdgeMaskVec)

    # Matrices for picking out ν and xy directions in derivatives 
    Pν  = spdiagm(vcat(ones(Int64, dimEdgeCount[1]), zeros(Int64, sum(dimEdgeCount[2:end]))))   # Diagonal sparse matrix to exclude all xy edges 
    Pxy  = spdiagm(vcat(zeros(Int64, dimEdgeCount[1]), ones(Int64, sum(dimEdgeCount[2:end]))))   # Diagonal sparse matrix to exclude all ν edges 

    # Weights
    W   = vertexVolumeWeightsMatrix(dimsPlus, spacing)
    W⁻¹ =  vertexVolumeWeightsInverseMatrix(dimsPlus, spacing)
    l⁻¹ = edgeLengthInverseMatrix(dimsPlus, spacing)

    ∇ₑ = l⁻¹*A       # Gradient operator giving gradient on each edge
    ∇cdot = -W⁻¹*Aᵀ  # Divergence operator giving divergence on each vertex calculated from edges 

    # Diffusivity field over edges 
    # Set no-flux boundary conditions by enforcing zero diffusivity in edges connection ghost points
    Aperpₑ = edgePerpendicularAreaMatrix(dimsPlus, spacing)
    𝓓ₑ     = 𝓓.*Aperpₑ # Sparse diagonal matrix of diffusivities over edges 

    # Diagonal matrices of compartment thickness h over all vertices hᵥ
    # Also diagonal matrix of thickness over edges, formed by taking mean of h at adjacent vertices 0.5.*Ā*hᵥ
    hᵥ_vec = reshape(mat_h, nVerts)         # Cisternal thickness evaluated over vertices 
    hₑ_vec = 0.5.*Ā*hᵥ_vec                  # Cisternal thickness evaluated over edges (mean of adjacent vertices)
    hᵥ = spdiagm(hᵥ_vec)                    # Cisternal thickness over vertices, as a sparse diagonal matrix
    hₑ = spdiagm(hₑ_vec)                    # Cisternal thickness over edges, as a sparse diagonal matrix
    aᵥ = spdiagm(1.0./(1.0 .+ α_C.*hᵥ_vec)) # Prefactor 1/(1+α_C*hᵥ(x)) evaluated over vertices, packaged into a sparse diagonal matrix for convenience
    aₑ = spdiagm(1.0./(1.0 .+ α_C.*hₑ_vec)) # Prefactor 1/(1+α_C*hₑ(x)) evaluated over edges, packaged into a sparse diagonal matrix for convenience

    u0fun(xs, μs, σs) = exp(-sum((xs.-μs).^2.0./σs.^2.0)) # Multidimensional Gaussian
    uMat = zeros(Float64, dimsPlus...)
    for ind in CartesianIndices(uMat)
        uMat[ind] = u0fun([νs[ind[1]]], [0.0], [νMax/10.0])
    end
    u0 = reshape(uMat, nVerts)
    # # For integration to normalise the number of monomers, we need to multiply the concentration at each point by the ν value of that point
    # νMat = ones(Int64, dimsPlus...)
    # for ii=2:dimsPlus[1]
    #     selectdim(νMat, 1, ii) .*= (ii-1)
    # end
    # νSparse = spdiagm(reshape(νMat, nVerts))
    # integ = sum(ghostVertexMaskSparse*W*νSparse*u0)
    # u0 .*= 𝓒/integ

    # Set value of Fₑ at each point in space
    matFₑTmp = ones(Float64, dimsPlus...)
    for i=2:length(dimsPlus)
        selectdim(matFₑTmp, i, 1) .*= 0.5
        selectdim(matFₑTmp, i, dimsPlus[i]) .*= 0.5
    end
    integF = prod(spacing[Not(1)])*sum(selectdim(matFₑTmp, 1, 1))
    # Ensure integral of Fₑ over space is π
    matFₑ = (π/integF).*ones(Float64, dimsPlus[Not(1)]...)
    matE = zeros(dimsPlus...)
    Esparse = spzeros(nVerts, nVerts)
    E!(u0, dimsPlus, Esparse, matE, matFₑ, K₂, dν)

    # PDE operator components
    L1 = aᵥ*∇cdot*Aperpₑ*(K₂*K₄.*Pν*∇ₑ - β.*Pν*Aᵤₚ)
    L2 = aᵥ*∇cdot*(hₑ*Pxy*𝓓ₑ*∇ₑ)
    p = (L1=L1, L2=L2, u0=u0, dimsPlus=dimsPlus, Esparse=Esparse, matE=matE, matFₑ=matFₑ, K₂=K₂, dν=dν)
    L = MatrixOperator(Esparse*L1.+L2, update_func! = updateOperator!)

    prob = ODEProblem(L, u0, (0.0, Tᵣ), p)
    sol = solve(prob, Vern9(), saveat=Tᵣ/100.0, progress=true)

    return sol
end

export glycosylationAnyD

end