
module UsefulFunctions

using LinearAlgebra
using SparseArrays


# function updateEintegral!(E, u, dimsPlus, K₂, matE, matFₑ)
#     cs = reshape(u, dimsPlus)     
#     for j = 1:dimsPlus[2]
#         integrationFactor = K₂/(K₂ + trapeziumRule(cs[:,j]))
#         matE[:,j] .= matFₑ[:,j].*integrationFactor
#     end
#     E .= spdiagm(reshape(matE, nVerts))
# end


function E!(u, dimsPlus, Esparse, matE, matFₑ, K₂, dν)
    # Convert state vector to matrix of concentrations (We're calculating enzyme distribution, but using bulk concentration?)
    cs = reshape(u, dimsPlus)
    # dν.*sum(cs[2:end-1,:], dims=1)[1,:] Gives integral of concentration over ν at each point in x
    matE .= matFₑ.*(K₂./(K₂ .+ dν.*sum(cs[2:end-1,:], dims=1)[1,:]))
    Esparse[diagind(Esparse)] .= repeat(matE, inner=dimsPlus[1])
    return nothing
end

# Integral of h*C over space 
function M_tilde(u, W, ghostVertexMaskVec, dims, dν, hᵥ)
    uInternal = reshape((W*hᵥ*u)[ghostVertexMaskVec], dims)
    M_tilde = sum(uInternal, dims=(2:length(dims)))
    return M_tilde./dν
end

# Dimensional bulk functional mass integrated over space and polymerisation 
function M_star(u, W, ghostVertexMaskVec, dims, hᵥ, ϕ, α_C, C_b, Ω)
    uInternal = reshape((W*hᵥ*u)[ghostVertexMaskVec], dims)
    M_tilde = sum(selectdim(uInternal, 1, round(Int, ϕ*dims[1]):dims[1]))
    prefactor = α_C*C_b*Ω/(π*(1+α_C))
    return prefactor*M_tilde
end

P_star(u, W, ghostVertexMaskVec, dims, hᵥ, ϕ, α_C, C_b, Ω, Tᵣstar) = M_star(u, W, ghostVertexMaskVec, dims, hᵥ, ϕ, α_C, C_b, Ω)/Tᵣstar

# function Mxy(u, W, ghostVertexMaskVec, dims, ϕ)
#     uInternal = reshape((W*u[ghostVertexMaskVec]), dims)
#     return sum(selectdim(uInternal, 1, round(Int, ϕ*dims[1]):end), dims=1)
# end

# function M(u, W, ghostVertexMaskVec, dims, ϕ)
#     uInternal = reshape((W*u[ghostVertexMaskVec]), dims)
#     return sum(selectdim(uInternal, 1, round(Int, ϕ*dims[1]):end))
# end

# function P(u, W, ghostVertexMaskVec, dims, ϕ, T)
#     uInternal = reshape((W*u[ghostVertexMaskSparse]), dims)
#     return sum(selectdim(uInternal, 1, round(Int, ϕ*dims[1]):end))/T
# end 


export M_tilde
export M_star
export P_star
export E!

end