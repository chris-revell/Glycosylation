#
#  Glycosylation.jl
#  Glycosylation
#

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
# ċ = a*(E*∇⋅(K₂*K₄.*Pν*∇ₑ - β*Pν*Aᵤₚ) + 𝒟.*∇⋅(hₑ*Pxy*∇ₑ))*cᵥ

# L = -W⁻¹*Aᵀ*𝒟*L⁻¹*A .+ W⁻¹*Aᵀ*V*Aᵤₚ # Express model as a matrix operator 


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
# ċ = a*(E*∇⋅(K₂*K₄.*Pν*∇ₑ - β*Pν*Aᵤₚ) + 𝒟.*∇⋅(hₑ*Pxy*∇ₑ))*cᵥ



module Glycosylation

using OrdinaryDiffEq
using SparseArrays
using UnPack
using FromFile
using DrWatson
using SciMLOperators
using Statistics
using InvertedIndices
using GaussianRandomFields
using LinearAlgebra

@from "$(srcdir("MakeIncidenceMatrix.jl"))" using MakeIncidenceMatrix
@from "$(srcdir("MakeWeightMatrices.jl"))" using MakeWeightMatrices
@from "$(srcdir("UsefulFunctions.jl"))" using UsefulFunctions
@from "$(srcdir("DerivedParameters.jl"))" using DerivedParameters

u0Fun(xs, μs, σs) = exp(-sum((xs.-μs).^2.0./σs.^2.0)) # Multidimensional Gaussian

function hFun(dims; λ=0.1, σ=0.1)
    if length(dims) == 2
        # cov = CovarianceFunction(length(dims)-1, Exponential(λ, σ=σ))# Gaussian(λ, σ=σ))
        cov = CovarianceFunction(length(dims)-1, Gaussian(λ, σ=σ))
        pts = range(0, stop=1, length=dims[2])
        grf = GaussianRandomField(cov, CirculantEmbedding(), pts, minpadding=113)
        # grf = GaussianRandomField(cov, Cholesky(), pts, minpadding=10001)
        mat_hSlice = sample(grf)[1:dims[2]]
    else
        # cov = CovarianceFunction(length(dims)-1, Exponential(λ, σ=σ))# Gaussian(λ, σ=σ))
        cov = CovarianceFunction(length(dims)-1, Gaussian(λ, σ=σ))
        pts1 = range(0, stop=sqrt(π), length=dims[2])
        pts2 = range(0, stop=sqrt(π), length=dims[3])
        grf = GaussianRandomField(cov, CirculantEmbedding(), pts1, pts2, minpadding=113)
        # grf = GaussianRandomField(cov, Cholesky(), pts1, pts2, minpadding=10001)
        mat_hSlice = sample(grf)[1:dims[2], 1:dims[3]]
    end
    mat_hSlice .= mat_hSlice.-mean(mat_hSlice).+1.0
    mat_h = zeros(dims...)
    for i=1:dims[1]
        selectdim(mat_h, 1, i) .= mat_hSlice
    end
    return mat_h
end

function hFunGaussian(dims; σ=0.5, μ=0.5)
    xMax = sqrt(π)
    xs   = collect(range(0.0, xMax, dims[2]))
    σx = xMax*σ
    μx = xMax*μ
    if length(dims) == 2
        mat_hSlice = [1.1 + 1.0*exp(-(x-μx)^2/σx^2) for x in xs]
        mat_h = zeros(dims...)
        for i=1:dims[1]
            selectdim(mat_h, 1, i) .= mat_hSlice
        end
        mat_h .= mat_h./mean(mat_h)
        return mat_h
    else
        mat_h = zeros(dims...)
        for i=1:dims[1]
            for j=1:dims[2]
                # selectdim(mat_h, 1, i) .= [0.1 + exp(-(x-μx)^2/σx^2 - (xs[i]-μx)^2/σx^2 ) for x in xs]
                mat_h[i, j, :] .= [1.1 + 1.0*exp(-(x-μx)^2/σx^2 - (xs[j]-μx)^2/σx^2 ) for x in xs]
            end
        end
        mat_h .= mat_h./mean(mat_h)
        return mat_h        
    end
end

function conditionProgressInterval(u, t, integrator)
    integrator.dt > t%integrator.p.interval
end

affectPrintProgress(integrator) = println("$(floor(Int64,100.0*integrator.t/integrator.p.tMax))")

function conditionHalfProduction(u, t, integrator)
    M̃ϕ(u, integrator.p.W, integrator.p.dims, integrator.p.dν, integrator.p.hᵥ, 0.5) > 0.5*π  
end

function conditionNuWall(u, t, integrator)
    M̃ϕ(u, integrator.p.W, integrator.p.dims, integrator.p.dν, integrator.p.hᵥ, 0.8) > 0.5*π  
end

# function conditionNuWall(u, t, integrator)
#     uInternal = reshape(u, integrator.p.dims...)
#     findmax(uInternal)[2][1] == integrator.p.dims[1]
#     # findmax(M̃(u, integrator.p.W, integrator.p.dims, integrator.p.dν, integrator.p.hᵥ))[2].I[1] > 0.9*integrator.p.dims[1] ? true : false
# end

affectTerminate!(integrator) = terminate!(integrator, ReturnCode.Success)    

# cbHalfProduction = DiscreteCallback(conditionHalfProduction, affectTerminate!)
cbHalfProduction = DiscreteCallback(conditionHalfProduction, affectTerminate!)
cbNuWall = DiscreteCallback(conditionNuWall, affectTerminate!)
cbProgress = DiscreteCallback(conditionProgressInterval, affectPrintProgress)

# Integrate over ν to find E field in spatial dimensions.
# When state vector u is reshaped to an array with shape dims, assume ν is the first dimension of this array
# Function is agnostic about the whether dims is of length 2 or 3.
function E!(u, dims, Esparse, matE, matFₑ, K₂, dν)
    uMat = reshape(u, dims...)
    integ = dν.*(sum(uMat, dims=1) .- 0.5.*selectdim(uMat, 1, 1) .- 0.5.*selectdim(uMat, 1, 1))
    # integ = dν.*(0.5.*selectdim(uMat, 1, 1) .+ dropdims(sum(selectdim(uMat, 1, 2:dims[1]-1), dims=1), dims=1) .+ 0.5.*selectdim(uMat, 1, dims[1]))
    for slice in eachslice(matE, dims=1)
        slice .= matFₑ.*(K₂./(K₂ .+ selectdim(integ, 1, 1)))
    end
    Esparse[diagind(Esparse)] .= reshape(matE, prod(dims))
    return nothing
end

# Function to update linear operator with new values for E at each iteration in solving the ODE system
function updateOperator!(L, u, p, t)
    # @unpack Part1, Part2, u0, dims, Esparse, matE, matFₑ, K₂, dν = p
    E!(u, p.dims, p.Esparse, p.matE, p.matFₑ, p.K₂, p.dν)
    L .= p.Esparse*p.Part1 .+ p.Part2
end

function glycosylation(dims, K₂, K₄, T̃ᵣ, α_C, 𝒟, β; thickness="uniform", fDist="uniform", differencing="centre", solver=SSPRK432(), nOutputs=100, λGRF=0.1, σGRF=0.1, σGaussian=0.1, μGaussian=0.5, terminateAt="T̃ᵣ", saveIntermediate=true)

    # PDE discretisation parameters 
    nSpatialDims = length(dims)-1
    
    xMax = sqrt(π)
    xs   = collect(range(0.0, xMax, dims[2]))
    dx   = xs[2]-xs[1]    
    yMax = xMax
    # ys   = collect(range(0.0, yMax, dims[3]))
    # dy   = ys[2]-ys[1]
    νMax = 1.0
    νs   = collect(range(0.0, νMax, dims[1])) # Positions of discretised vertices in polymerisation space 
    dν   = νs[2]-νs[1]
    nSpatialDims == 1 ? spacing  = [dν, dx, xMax] : spacing  = [dν, dx, dx]

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
    𝒟ₑ     = 𝒟.*Aperpₑ # Sparse diagonal matrix of diffusivities over edges 

    # Diagonal matrices of compartment thickness h over all vertices hᵥ
    # Also diagonal matrix of thickness over edges, formed by taking mean of h at adjacent vertices 0.5.*Ā*hᵥ
    if thickness=="GRF"
        mat_h = hFun(dims, λ=λGRF, σ=σGRF)
    elseif thickness=="Gaussian"
        mat_h = hFunGaussian(dims, σ=σGaussian, μ=μGaussian)
    else 
        mat_h = ones(dims...)
    end
    hᵥ_vec = reshape(mat_h, nVerts)         # Cisternal thickness evaluated over vertices 
    hₑ_vec = 0.5.*Ā*hᵥ_vec                  # Cisternal thickness evaluated over edges (mean of adjacent vertices)
    hᵥ = spdiagm(hᵥ_vec)                    # Cisternal thickness over vertices, as a sparse diagonal matrix
    hₑ = spdiagm(hₑ_vec)                    # Cisternal thickness over edges, as a sparse diagonal matrix
    aᵥ = spdiagm(1.0./(1.0 .+ α_C.*hᵥ_vec)) # Prefactor 1/(1+α_C*hᵥ(x)) evaluated over vertices, packaged into a sparse diagonal matrix for convenience

    uMat = zeros(Float64, dims...)
    for ind in CartesianIndices(uMat)
        uMat[ind] = u0Fun([νs[ind[1]]], [0.0], [νMax/50.0])
    end
    # Ensure that the integral of concentration over ν at each point in space is 1
    integ = spacing[1].*(0.5.*selectdim(uMat, 1, 1) .+ dropdims(sum(selectdim(uMat, 1, 2:dims[1]-1), dims=1), dims=1) .+ 0.5.*selectdim(uMat, 1, dims[1]))    
    
    u0 = reshape(uMat, nVerts)
    u0 .*= 1.0/integ[1]
    
    # Set value of Fₑ at each point in space
    # Integral of Fₑ over space is π
    if fDist == "uniform"
        matFₑ = ones(Float64, dims[Not(1)]...)
        matFₑTmp = copy(matFₑ)
    else
        matFₑ = selectdim(hFunGaussian(dims, σ=σGaussian, μ=μGaussian), 1, 1)
        matFₑTmp = copy(matFₑ)
    end
    for i=1:length(size(matFₑTmp))
        selectdim(matFₑTmp, i, 1) .*= 0.5
        selectdim(matFₑTmp, i, size(matFₑTmp)[i]) .*= 0.5
    end
    integF = prod(spacing[Not(1)])*sum(selectdim(matFₑTmp, 1, 1:size(matFₑTmp)[1]))
    
    # Ensure integral of Fₑ over space is π
    # matFₑ = (1/integF).*ones(Float64, dims[Not(1)]...)
    matFₑ .*= π/integF
    matE = zeros(dims...)

    Esparse = spzeros(nVerts, nVerts)
    E!(u0, dims, Esparse, matE, matFₑ, K₂, spacing[1])

    if differencing=="upstream"
        Part1 = aᵥ*∇cdot*Aperpₑ*(K₂*K₄.*Pν*∇ₑ - β.*Pν*Aᵤₚ)
    else
        Part1 = aᵥ*∇cdot*Aperpₑ*(K₂*K₄.*Pν*∇ₑ - β.*Pν*Ā./2.0)
    end
    Part2 = aᵥ*∇cdot*Aperpₑ*(hₑ*Pxy*𝒟ₑ*∇ₑ)

    p = (Part1 = Part1, 
        Part2 = Part2, 
        u0 = u0, 
        dims = dims, 
        Esparse = Esparse, 
        matE = matE, 
        matFₑ = matFₑ, 
        K₂ = K₂, 
        dν = dν,
        W = W,
        hᵥ = hᵥ,
        interval = T̃ᵣ/(nOutputs-1),
        tMax = T̃ᵣ
    )
    fullOperator = MatrixOperator(Esparse*Part1, update_func! = updateOperator!)
    println("solving")
    if terminateAt == "halfProduction"
        prob = ODEProblem(fullOperator, u0, (0.0, T̃ᵣ), p)
        if saveIntermediate == true
            # sol = solve(prob, solver, tstops= callback=cbHalfProduction, save_on=false, save_start=false, save_end=true)#, dt=0.0001) , saveat=T̃ᵣ/(nOutputs-1)
            sol = solve(prob, solver, tstops=T̃ᵣ/(nOutputs-1), callback=cbHalfProduction, saveat=T̃ᵣ/(nOutputs-1), save_end=true) 
        else
            sol = solve(prob, solver, tstops=T̃ᵣ/(nOutputs-1), callback=cbHalfProduction, save_on=false, save_end=true) 
        end
    elseif terminateAt == "nuWall"
        prob = ODEProblem(fullOperator, u0, (0.0, T̃ᵣ), p)
        if saveIntermediate == true
            sol = solve(prob, solver, tstops=T̃ᵣ/(nOutputs-1), callback=CallbackSet(cbNuWall, cbProgress), saveat=T̃ᵣ/(nOutputs-1), save_end=true)
        else
            sol = solve(prob, solver, tstops=T̃ᵣ/(nOutputs-1), callback=CallbackSet(cbNuWall, cbProgress), save_on=false, save_end=true) 
        end
    else 
        prob = ODEProblem(fullOperator, u0, (0.0, T̃ᵣ), p)
        if saveIntermediate == true
            sol = solve(prob, solver, tstops=T̃ᵣ/(nOutputs-1), saveat=T̃ᵣ/(nOutputs-1), save_end=true)
        else
            sol = solve(prob, solver, tstops=T̃ᵣ/(nOutputs-1), save_on=false, save_end=true) 
        end
    end

    return sol, p
end

export glycosylation, hFun, hFunGaussian

end