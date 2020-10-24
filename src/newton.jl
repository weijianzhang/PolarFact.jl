#
# Newton's method
#
# Reference:
# Nicholas J. Higham, Computing the Polar Decomposition ---with Applications,
# SIAM J. Sci. Statist. Comput. Vol. 7, Num 4 (1986) pp. 1160-1174.
#
   
mutable struct NewtonAlg{T} <: PolarAlg
    maxiter::Int        # maximum number of iterations.
    scale::Bool         # whether to scale Newton iteration.
    verbose::Bool       # whether to show procedural information
    tol::T        # tolerance for convergence
    scale_tol::T  # tolerance for acceleration scaling 

    function NewtonAlg{T}( ;maxiter::Integer=100,
                       verbose::Bool=false,
                       scale::Bool=true,
                       tol::Real=cbrt(eps(T)),
                       scale_tol::Real=eps(T)^(1/4)) where T
        maxiter > 1 || error("maxiter must be greater than 1.")
        tol > 0 || error("tol must be positive.")
        scale_tol > 0 || error("scale_tol must be positive.")
        
        new{T}(maxiter,
            scale,
            verbose,
            tol,
            scale_tol)
    end
end


function solve!(alg::NewtonAlg{T}, 
                   X::Matrix{T}, U::Matrix{T}, H::Matrix{T}) where {T}
    common_iter_scal!(NewtonUpdater(alg.scale, alg.scale_tol), X, U, H, alg.maxiter, alg.verbose, alg.tol)
end

mutable struct NewtonUpdater{T} <: PolarUpdater 
    scale::Bool
    scale_tol::T
end

function update_U!(upd::NewtonUpdater, U::Matrix{T}) where {T}
    scale = upd.scale
    Uinv = Array{T}(undef, size(U))
    Uinvt = Array{T}(undef, size(U))
    copyto!(Uinv, inv(U))
    
    # 1, Inf-norm scaling 
    if scale
        g = (opnorm(Uinv,1) * opnorm(Uinv, Inf) / (opnorm(U,1) * opnorm(U, Inf)) )^(1/4)
    else
        g = one(T)
    end
    transpose!(Uinvt, Uinv)
    copyto!(U, (g * U + Uinvt /g) / convert(T, 2))

end



# Newton Schulz algrithm

mutable struct NewtonSchulzAlg{T} <: PolarAlg  
    maxiter::Int     # maximum number of iterations.
    verbose::Bool    # whether to show procedural information
    tol::T     # convergence tolerance. 

    function NewtonSchulzAlg{T}(; maxiter::Integer=100,
                             verbose::Bool=false,
                             tol::Real=cbrt(eps(T))) where T
        maxiter > 1 || error("maxiter must be greater than 1.")
        tol > 0 || error("tol must be positive.")

        new{T}(int(maxiter),
            verbose,
            tol)
    end
end

function solve!(alg::NewtonSchulzAlg, 
                   X::Matrix{T}, U::Matrix{T}, H::Matrix{T}) where {T}
    
    # Newton Schulz converge quadratically if norm(X) < sqrt(3)
    opnorm(X) < convert(T, sqrt(3)) || throw(ArgumentError("The norm of the input matrix must be smaller than sqrt(3)."))

    common_iter!(NewtonSchulzUpdater(), X, U, H, alg.maxiter, alg.verbose, alg.tol)
end

struct NewtonSchulzUpdater <: PolarUpdater end

function update_U!(upd::NewtonSchulzUpdater, U::Matrix{T}) where {T}
    UtU = Array{T}(undef, size(U))
    mul!(UtU, transpose(U), U)
    copyto!(U, 0.5*U*(3*I - UtU))
end
