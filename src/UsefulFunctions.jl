
module UsefulFunctions

using LinearAlgebra
using SparseArrays

simpsonsRule(xs) = sum(xs[1:end-1].+xs[2:end])./2.0

function productionTotalM(u, W, ghostVertexMask, dims, ϕ)
    uInternal = reshape((W*u)[ghostVertexMask], dims)
    return sum(selectdim(uInternal, 1, round(Int, ϕ*dims[1])))
end

function updateEintegral!(E, u, dimsPlus, K₂, matE, matFₑ)
    cs = reshape(u, dimsPlus)     
    for j = 1:dimsPlus[2]
        integrationFactor = K₂/(K₂ + simpsonsRule(cs[:,j]))
        matE[:,j] .= matFₑ[:,j].*integrationFactor
    end
    E .= spdiagm(reshape(matE, nVerts))
end

export simpsonsRule
export productionTotalM

end