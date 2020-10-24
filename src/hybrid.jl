#
# Newton and Newton Schulz Hybrid algoritm
#
# Reference:
# Nicholas J. Higham and Robert S. Schreiber, Fast Polar Decomposition
# of an arbitrary matrix, SIAM, J. Sci. Statist. Comput. Vol. 11, No. 4
# (1990) pp. 648-655.
#

mutable struct NewtonHybridAlg{T} <: PolarAlg
    maxiter::Int    # maximum number of iterations.
    verbose::Bool   # whether to show procedural information. 
    tol::T    # convergence tolerance
    theta::T  # switch parameter
    
    function NewtonHybridAlg{T}( ; maxiter::Integer=100,
                             verbose::Bool=false,
                             tol::Real=cbrt(eps(T)),
                             theta::Real=convert(T, 0.6)) where T
        maxiter > 1 || error("maxiter must  be greater than 1.")
        tol > 0 || error("tol must be positive.")
        theta > 0 || error("theta must be positive.")
        
        new{T}(maxiter, 
            verbose,
            tol,
            theta)
    end 
end
    
function solve!(alg::NewtonHybridAlg, 
                X::Matrix{T}, U::Matrix{T}, H::Matrix{T}) where {T}
    common_iter_hybr!(NewtonUpdater(true, eps(T)^(1/4)), NewtonSchulzUpdater(), X, U, H, 
                      alg.maxiter, alg.verbose, alg.tol, alg.theta)
end

