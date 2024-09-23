#
#  MakeWeightMatrices.jl
# Glycosylation
#
#  Created by Christopher Revell on 04/06/2024.

module MakeWeightMatrices

using SparseArrays
using LinearAlgebra
using InvertedIndices

# function solChecks(sol, W, ghostVertexMaskSparse)
#     mass1 = sum((W*sol.u[1])[ghostVertexMaskSparse])
#     massEnd = sum((W*sol.u[end])[ghostVertexMaskSparse])
#     @show mass1
#     @show massEnd
#     minima = [minimum(u[ghostVertexMaskSparse]) for u in sol.u]
#     all_t_min = minimum(minima)
#     @show all_t_min
#     return nothing
# end

function solChecks(sol, W, ghostVertexMaskSparse)
    mass1 = sum(ghostVertexMaskSparse*W*sol.u[1])
    massEnd = sum(ghostVertexMaskSparse*W*sol.u[end])
    @show mass1
    @show massEnd
    minima = [minimum(ghostVertexMaskSparse*u) for u in sol.u]
    all_t_min = minimum(minima)
    @show all_t_min
    return nothing
end


# Matrices for picking out ν and xy directions in derivatives 
# Matrix of i-directed edge accessibility
# P_i = ones(Nxplus-1, Nνplus); P_i[:, 1] .= 0.0; P_i[:, end] .= 0.0; P_i[1, :] .= 0.0; P_i[end, :] .= 0.0
# # Matrix of j-directed edge accessibility  
# P_j = ones(Nxplus, Nνplus-1); P_j[:, 1] .= 0.0; P_j[:, end] .= 0.0; P_j[1, :] .= 0.0; P_j[end, :] .= 0.0
# # Matrix of k-directed edge accessibility  
# P_j = ones(Nxplus, Nνplus-1); P_j[:, 1] .= 0.0; P_j[:, end] .= 0.0; P_j[1, :] .= 0.0; P_j[end, :] .= 0.0
# P = dropzeros(spdiagm(vcat(reshape(P_i, nEdgesi), reshape(P_j, nEdgesj)))) # Diagonal sparse matrix to exclude all edges adjacent to ghost points  
# Pν = dropzeros(spdiagm(vcat(zeros(nEdgesi), reshape(P_j, nEdgesj)))) # Diagonal sparse matrix to exclude all xy edges and ν edges adjacent to ghost points  
# Pxy = dropzeros(spdiagm(vcat(reshape(P_i, nEdgesi), zeros(nEdgesj)))) # Diagonal sparse matrix to exclude all ν edges and xy edges adjacent to ghost points 


# Ghost point mask is a 1D vector in which the value at component i 
# is true if vertex i in the flattened vector of vertices is an internal vertex
# but false if vertex i is a ghost vertex in the flattened vector of vertices 
# This can be used to exclude ghost points from calculations over the whole state vector
function makeGhostVertexMask(dimsPlus)
    ghostMaskVertex = fill(true, dimsPlus...)
    for i=1:length(dimsPlus)
        selectdim(ghostMaskVertex, i, 1) .= false
        selectdim(ghostMaskVertex, i, dimsPlus[i]) .= false
    end
    # ghostMaskVertex = spdiagm(reshape(ghostMaskVertex, prod(dimsPlus)))
    return reshape(ghostMaskVertex, prod(dimsPlus))
end

function makeGhostEdgeMask(dimsPlus)
    ghostMaskEdge = Bool[]
    for i=1:length(dimsPlus)
        dimsVec_i = copy(dimsPlus)
        dimsVec_i[i] = dimsVec_i[i]-1
        nEdgesi = prod(dimsVec_i)
        l_i = fill(true, dimsVec_i...)
        l_iSize = size(l_i)
        for (j,s) in enumerate(size(l_i))
            selectdim(l_i, j, 1) .= false
            selectdim(l_i, j, s) .= false
        end
        append!(ghostMaskEdge, reshape(l_i, nEdgesi))
    end
    return ghostMaskEdge
end


# function makeGhostEdgeMaskNew(dimsPlus)
#     dimEdgeCount = Int64[]
#     for i=1:length(dimsPlus)
#         push!(dimEdgeCount, (dimsPlus[i]-1)*prod(dimsPlus[Not(i)]))
#     end
#     nEdges  = sum(dimEdgeCount)

#     ghostEdgeMaskVec = fill(true, nEdges)

#     dimEdgeArrayStrides = (1, (dimsPlus[1]-1), (dimsPlus[1]-1)*dimsPlus[2])
#     for kk=1:dimsPlus[3]
#         for jj=1:dimsPlus[2]
#             for ii in [1,(dimsPlus[1]-1)]
#                 edgeIndex = 1 + (ii-1)*dimEdgeArrayStrides[1] + (jj-1)*dimEdgeArrayStrides[2] + (kk-1)*dimEdgeArrayStrides[3]
#                 ghostEdgeMaskVec[edgeIndex] = false
#             end
#         end
#     end  
#     dimEdgeArrayStrides = (1, dimsPlus[1], dimsPlus[1]*(dimsPlus[2]-1))
#     for kk=1:dimsPlus[3]
#         for jj in [1,(dimsPlus[2]-1)]
#             for ii=1:dimsPlus[1]
#                 edgeIndex = dimEdgeCount[1] + 1 + (ii-1)*dimEdgeArrayStrides[1] + (jj-1)*dimEdgeArrayStrides[2] + (kk-1)*dimEdgeArrayStrides[3]
#                 ghostEdgeMaskVec[edgeIndex] = false
#             end
#         end
#     end  
#     dimEdgeArrayStrides = (1, dimsPlus[1], dimsPlus[1]*dimsPlus[2])
#     for kk in [1,(dimsPlus[3]-1)]
#         for jj=1:dimsPlus[2]
#             for ii=1:dimsPlus[1]
#                 edgeIndex = dimEdgeCount[1] + dimEdgeCount[2] + 1 + (ii-1)*dimEdgeArrayStrides[1] + (jj-1)*dimEdgeArrayStrides[2] + (kk-1)*dimEdgeArrayStrides[3]
#                 ghostEdgeMaskVec[edgeIndex] = false
#             end
#         end
#     end  
#     return ghostEdgeMaskVec
# end




# Vertex weights
# Forming a diagonal matrix of volumes around each vertex, 
# divided by 2 at the periphery
function vertexVolumeWeightsMatrix(dims, spacing)
    matW = fill(prod(spacing), dims...)
    for i=1:length(dims)
        selectdim(matW, i, 1) ./= 2.0
        selectdim(matW, i, dims[i]) ./= 2.0
    end
    vecW = reshape(matW, prod(dims))
    W = spdiagm(vecW)
    return W
end

# Vertex weights inverse
# Forming a diagonal matrix of 1/volumes around each vertex, 
# with volumes divided by 2 at the periphery
function vertexVolumeWeightsInverseMatrix(dims, spacing)
    matW = fill(prod(spacing), dims...)
    for i=1:length(dims)
        selectdim(matW, i, 1) ./= 2.0
        selectdim(matW, i, dims[i]) ./= 2.0
    end
    vecW = reshape(matW, prod(dims))
    W⁻¹ = spdiagm(1.0./vecW)
    return W⁻¹
end

# Diagonal matrix of edge lengths
function edgeLengthMatrix(dims, spacing)
    lvec = Float64[]
    for i=1:length(dims)
        nEdgesi = (dims[i]-1)*prod(dims[Not(i)])
        l_i = fill(spacing[i], nEdgesi)
        append!(lvec, l_i)
    end
    l = spdiagm(lvec)
    return l
end

# Diagonal inverse matrix of edge lengths
function edgeLengthInverseMatrix(dims, spacing)
    lvec = Float64[]
    for i=1:length(dims)
        nEdgesi = (dims[i]-1)*prod(dims[Not(i)])
        l_i = fill(spacing[i], nEdgesi)
        append!(lvec, l_i)
    end
    l⁻¹ = spdiagm(1.0./lvec)
    return l⁻¹
end

# Diagonal matrix of areas perpendicular to each edge, 
# meaning the area through which diffusive 
# flux in the direction of a given edge passes
# Factor of 1/2 applied to peripheral edges, assuming peripheral 
# vertices lie on the edge of the solution domain so that the surrounding
# grid points are halved.
function edgePerpendicularAreaMatrix(dims, spacing)  
    Aperpvec = Float64[]
    for i=1:length(dims)
        edgeDims = copy(dims)
        edgeDims[i] -= 1
        nEdgesi = prod(edgeDims)
        AMat = fill(prod(spacing[Not(i)]), edgeDims...)
        for j = [jj for jj in 1:length(edgeDims) if jj!=i]
            selectdim(AMat, j, 1) ./= 2.0
            selectdim(AMat, j, edgeDims[j]) ./= 2.0
        end
        append!(Aperpvec, reshape(AMat, nEdgesi))
    end
    Aperpₑ   = spdiagm(Aperpvec) 
    return Aperpₑ
end

function setBoundaryEdgesToZero(dims)
    ghostMaskArray = fill(true, dims...)
    for i=1:length(dims)
        selectdim(ghostMaskArray1, i, 1) .= false
        selectdim(ghostMaskArray1, i, dims[i]) .= false
    end
    ghostMask = reshape(ghostMaskArray, prod(dims))
    return ghostMask
end

export solChecks
export makeGhostVertexMask
export makeGhostEdgeMask
export makeGhostEdgeMaskNew
export vertexVolumeWeightsMatrix
export vertexVolumeWeightsInverseMatrix
export edgeLengthMatrix
export edgeLengthInverseMatrix
export edgePerpendicularAreaMatrix
export setBoundaryEdgesToZero

end

# # Diagonal matrix of volumes around each edge, divided by 2 at the periphery
# # Matrix of i-directed edge weights  
# F_i = fill(dx*dy*dν, (Nxplus-1, Nyplus, Nνplus))
# F_i[:, :, 1] ./= 2.0
# F_i[:, :, end] ./= 2.0
# F_i[:, 1, :] ./= 2.0
# F_i[:, end, :] ./= 2.0
# F_i[1, :, :] ./= 2.0
# F_i[end, :, :] ./= 2.0
# # Matrix of j-directed edge weights  
# F_j = fill(dx*dy*dν, (Nxplus, Nyplus-1, Nνplus))
# F_j[:, :, 1] ./= 2.0
# F_j[:, :, end] ./= 2.0
# F_j[:, 1, :] ./= 2.0
# F_j[:, end, :] ./= 2.0
# F_j[1, :, :] ./= 2.0
# F_j[end, :, :] ./= 2.0
# # Matrix of k-directed edge weights  
# F_k = fill(dx*dy*dν, (Nxplus, Nyplus, Nνplus-1))
# F_k[:, :, 1] ./= 2.0
# F_k[:, :, end] ./= 2.0
# F_k[:, 1, :] ./= 2.0
# F_k[:, end, :] ./= 2.0
# F_k[1, :, :] ./= 2.0
# F_k[end, :, :] ./= 2.0
# Fvec = vcat(reshape(F_i, nEdgesi), reshape(F_j, nEdgesj), reshape(F_k, nEdgesk))
# F = spdiagm(Fvec)
# F⁻¹ = spdiagm(1.0./Fvec)

# For creating velocity field, diffusivity field, etc with peripheral edges set to zero
# function dimensionSpecificEdgeWeights
#     # Velocity field 
#     V_i = fill(0.0, (Nxplus-1, Nyplus, Nνplus))
#     V_j = fill(0.0, (Nxplus, Nyplus-1, Nνplus))
#     V_k = fill(0.0, (Nxplus, Nyplus, Nνplus-1))
#     for k=2:Nνplus-2
#         for j=2:Nyplus-1
#             for i=2:Nxplus-1
#                 V_k[i,j,k] = β/(1+α*hFun(xs[i], 0.5, 0.1, ys[j], 0.5, 0.1))
#             end
#         end
#     end
#     Vvec = vcat(reshape(V_i, nEdgesi), reshape(V_j, nEdgesj), reshape(V_k, nEdgesk))
#     V = spdiagm(Vvec)   # Diagonal matrix of advection velocities at each edge

# # Diffusivity field over edges 
# # Set no-flux boundary conditions by enforcing zero diffusivity in edges connection ghost points
# D_i = fill(0.0, (Nxplus-1, Nyplus, Nνplus))
# D_j = fill(0.0, (Nxplus, Nyplus-1, Nνplus))
# D_k = fill(0.0, (Nxplus, Nyplus, Nνplus-1))
# for k=2:Nνplus-2
#     for j=2:Nyplus-1
#         for i=2:Nxplus-1
#             D_k[i,j,k] = 0.11 #dx*dy*K₂*K₄/(1+α*hFun(xs[i], 0.5, 0.1, ys[j], 0.5, 0.1))
#         end
#     end
# end
# for k=2:Nνplus-1
#     for j=2:Nyplus-2
#         for i=2:Nxplus-1
#             D_j[i,j,k] = 0.12 #dx*dν*K₂*K₄/(1+α*hFun(xs[i], 0.5, 0.1, ys[j], 0.5, 0.1))
#         end
#     end
# end
# for k=2:Nνplus-1
#     for j=2:Nyplus-1
#         for i=2:Nxplus-2
#             D_i[i,j,k] = 0.13 #dy*dν*K₂*K₄/(1+α*hFun(xs[i], 0.5, 0.1, ys[j], 0.5, 0.1))
#         end
#     end
# end
# Dvec = vcat(reshape(D_i, nEdgesi), reshape(D_j, nEdgesj), reshape(D_k, nEdgesk))
# D = spdiagm(Dvec) # Diagonal matrix of advection velocities at each edge