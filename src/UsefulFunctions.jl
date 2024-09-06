
module UsefulFunctions

using LinearAlgebra
using SparseArrays
using UnPack

# Integrate over ν to find E field in spatial dimensions.
# When state vector u is reshaped to an array with shape dimsPlus, assume ν is the first dimension of this array
# Function is agnostic about the whether dimsPlus is of length 2 or 3.
function E!(u, dimsPlus, Esparse, matE, matFₑ, K₂, dν)
    # Convert state vector to matrix of concentrations (We're calculating enzyme distribution, but using bulk concentration?)
    cs = selectdim(reshape(u, dimsPlus...), 1, 2:(dimsPlus[1]-1))
    # dν.*sum(cs[2:end-1,:], dims=1)[1,:] Gives integral of concentration over ν at each point in x
    matE .= matFₑ.*(K₂./(K₂ .+ dν.*dropdims(sum(cs, dims=1), dims=1)))
    Esparse[diagind(Esparse)] .= repeat(reshape(matE, prod(size(matE))), inner=dimsPlus[1])
    return nothing
end

# Function to update linear operator with new values for E at each iteration in solving the ODE system
function updateOperator!(L, u, p, t)
    @unpack L1, L2, u0, dimsPlus, Esparse, matE, matFₑ, K₂, dν = p
    E!(u, dimsPlus, Esparse, matE, matFₑ, K₂, dν)
    L .= Esparse*L1 .+ L2
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

export M_tilde
export M_star
export P_star
export E!
export updateOperator!

end