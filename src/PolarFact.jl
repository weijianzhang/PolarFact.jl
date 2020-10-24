module PolarFact

using LinearAlgebra
using Printf

import LinearAlgebra: transpose!

export polarfact


include("common.jl")     # common functions and types
include("newton.jl")     # using Newton's method
include("svd.jl")        # using the SVD method
include("halley.jl")     # using Halley's method and the QDWH method
include("hybrid.jl")     # using hybrid method

function polarfact(A::Matrix{T}; 
                   alg::Symbol=:newton, 
                   maxiter::Integer=:100,
                   tol::Real = cbrt(eps(T)),
                   verbose::Bool=false) where {T}

    # choose algorithm 
    algorithm = 
       alg == :newton ? NewtonAlg{T}(; maxiter=maxiter, tol=tol, verbose=verbose) :
       alg == :qdwh ?  QDWHAlg{T}(; maxiter=maxiter, tol=tol, verbose=verbose) : 
       alg == :halley ? HalleyAlg{T}(; maxiter=maxiter, tol=tol, verbose=verbose) :
       alg == :svd ? SVDAlg() :   
       alg == :schulz ? NewtonSchulzAlg{T}(; maxiter=maxiter, tol=tol, verbose=verbose) :
       alg == :hybrid ? NewtonHybridAlg{T}(; maxiter=maxiter, tol=tol, verbose=verbose) :
       error("Invalid algorithm.")

    # Initialization: if m > n, do QR factorization 
    m, n = size(A)
    if m > n
        if alg in [:newton, :qdwh, :halley, :schulz, :hybrid]
            m = n
            A = qrfact(A)[:R]
        end
    elseif m < n
        error("The row dimension of the input matrix must be 
              greater or equal to column dimension.")
    end
    U = Array{T}(undef, m, n)
    H = Array{T}(undef, n, n)
    # solve for polar factors
    solve!(algorithm, A, U, H)
end

end # module
