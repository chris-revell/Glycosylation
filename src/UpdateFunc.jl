#
#  UpdateFunc.jl
#  GolgiModels
#
#  Created by Christopher Revell on 04/06/2024.

module UpdateFunc

using SparseArrays
using LinearAlgebra
using InvertedIndices

function update_func!(L, u, p, t)
    @unpack L1,
        L2,
        Nνplus,
        Nxplus,
        Nyplus,
        K₂,
        matE,
        E,
        matFₑ = p

    cs = reshape(u, (Nνplus, Nxplus))     
    for j = 1:Nxplus
        integrationFactor = K₂/(K₂ + simpsonsRule(cs[:,j]))
        matE[:,j] .= matFₑ[:,j].*integrationFactor
    end
    E .= spdiagm(reshape(matE, nVerts)) 
    L .= E*L1 .+ L2
end