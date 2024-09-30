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

# L = -W⁻¹*Aᵀ*𝓓*L⁻¹*A .+ W⁻¹*Aᵀ*V*Aᵤₚ # Express model as a matrix operator 


# Cνν = W⁻¹*Aᵀ*Pν*L⁻¹*A
# Cν = Aᵀ*L⁻¹*Pν*Aᵤₚ
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
@from "$(srcdir("UsefulFunctions.jl"))" using UsefulFunctions
@from "$(srcdir("DerivedParameters.jl"))" using DerivedParameters


function glycosylationAnyD(mat_h, dims, Ωperp, 𝓒, K₂, K₄, Tᵣ, α_C , 𝓓, β)

    # PDE discretisation parameters 
    nSpatialDims = length(dims)-1
    
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

    h₀ = mean(selectdim(mat_h, 1, 1))

    A   = makeIncidenceMatrix3D(dims)
    Ā   = abs.(A)
    Aᵀ  = transpose(A)
    Aᵤₚ = dropzeros((Ā-A).÷2)

    # Number of edges over each dimension 
    dimEdgeCount = Int64[]
    for i=1:length(dims)
        push!(dimEdgeCount, (dims[i]-1)*prod(dims[Not(i)]))
    end
    nVerts  = prod(dims)          # Total number of vertices 
    nEdges  = sum(dimEdgeCount)   # Total number of edges over all dimensions 

    # Matrices for picking out ν and xy directions in derivatives 
    Pν  = spdiagm(vcat(ones(Int64, dimEdgeCount[1]), zeros(Int64, sum(dimEdgeCount[2:end]))))   # Diagonal sparse matrix to exclude all xy edges 
    Pxy  = spdiagm(vcat(zeros(Int64, dimEdgeCount[1]), ones(Int64, sum(dimEdgeCount[2:end]))))   # Diagonal sparse matrix to exclude all ν edges 

    # Weights
    W   = vertexVolumeWeightsMatrix(dims, spacing)
    W⁻¹ =  vertexVolumeWeightsInverseMatrix(dims, spacing)
    L⁻¹ = edgeLengthInverseMatrix(dims, spacing)

    ∇ₑ = L⁻¹*A       # Gradient operator giving gradient on each edge
    ∇cdot = -W⁻¹*Aᵀ  # Divergence operator giving divergence on each vertex calculated from edges 

    # Diffusivity field over edges 
    # Set no-flux boundary conditions by enforcing zero diffusivity in edges connection ghost points
    Aperpₑ = edgePerpendicularAreaMatrix(dims, spacing)
    𝓓ₑ     = 𝓓.*Aperpₑ # Sparse diagonal matrix of diffusivities over edges 

    # Diagonal matrices of compartment thickness h over all vertices hᵥ
    # Also diagonal matrix of thickness over edges, formed by taking mean of h at adjacent vertices 0.5.*Ā*hᵥ
    hᵥ_vec = reshape(mat_h, nVerts)         # Cisternal thickness evaluated over vertices 
    hₑ_vec = 0.5.*Ā*hᵥ_vec                  # Cisternal thickness evaluated over edges (mean of adjacent vertices)
    hᵥ = spdiagm(hᵥ_vec)                    # Cisternal thickness over vertices, as a sparse diagonal matrix
    hₑ = spdiagm(hₑ_vec)                    # Cisternal thickness over edges, as a sparse diagonal matrix
    aᵥ = spdiagm(1.0./(1.0 .+ α_C.*hᵥ_vec)) # Prefactor 1/(1+α_C*hᵥ(x)) evaluated over vertices, packaged into a sparse diagonal matrix for convenience

    u0fun(xs, μs, σs) = exp(-sum((xs.-μs).^2.0./σs.^2.0)) # Multidimensional Gaussian
    uMat = zeros(Float64, dims...)
    for ind in CartesianIndices(uMat)
        uMat[ind] = u0fun([νs[ind[1]]], [0.0], [νMax/10.0])
    end
    u0 = reshape(uMat, nVerts)
    integ = sum(W*u0)
    u0 .*= 𝓒/integ
    
    # Set value of Fₑ at each point in space
    matFₑTmp = ones(Float64, dims...)
    for i=2:length(dims)
        selectdim(matFₑTmp, i, 1) .*= 0.5
        selectdim(matFₑTmp, i, dims[i]) .*= 0.5
    end
    integF = prod(spacing[Not(1)])*sum(selectdim(matFₑTmp, 1, 1))
    # Ensure integral of Fₑ over space is π
    matFₑ = (π/integF).*ones(Float64, dims[Not(1)]...)
    matE = zeros(dims...)
    Esparse = spzeros(nVerts, nVerts)
    E!(u0, dims, Esparse, matE, matFₑ, K₂, dν)

    # PDE operator components
    # Part1 = aᵥ*∇cdot*(K₂*K₄.*Pν*∇ₑ - β.*Pν*Aᵤₚ)
    # Part2 = 𝓓.*aᵥ*∇cdot*(hₑ*Pxy*∇ₑ)

    Part1 = aᵥ*∇cdot*Aperpₑ*(K₂*K₄.*Pν*∇ₑ - β.*Pν*Aᵤₚ)
    Part2 = aᵥ*∇cdot*(hₑ*Pxy*𝓓ₑ*∇ₑ)

    p = (Part1 = Part1, 
        Part2 = Part2, 
        u0 = u0, 
        dims = dims, 
        Esparse = Esparse, 
        matE = matE, 
        matFₑ = matFₑ, 
        K₂ = K₂, 
        dν = dν,
    )
    fullOperator = MatrixOperator(Esparse*Part1.+Part2, update_func! = updateOperator!)
    prob = ODEProblem(fullOperator, u0, (0.0, Tᵣ), p)
    println("solving")
    sol = solve(prob, Vern9(), saveat=Tᵣ/50.0, progress=true)

    return sol
end

export glycosylationAnyD

end