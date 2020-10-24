# 
# Halley's method
# 
# Reference:
# Y. Nakatsukasa, Z. Bai and F. Gygi, Optimizing Halley's iteration 
# for computing the matrix polar decomposition, SIAM, J. Mat. Anal. 
# Vol. 31, Num 5 (2010) pp. 2700-2720 
#

mutable struct HalleyAlg{T} <: PolarAlg
    maxiter::Int
    verbose::Bool
    tol::T
    
    function HalleyAlg{T}( ;maxiter::Integer=100,
                       verbose::Bool=false,
                       tol::Real = cbrt(eps(T))) where T
        maxiter > 1 || error("maxiter must be greater than 1.")
        tol > 0 || error("tol must be positive.")
        
        new{T}(maxiter,
            verbose,
            tol)
    end
end

function solve!(alg::HalleyAlg{T},
                X::Matrix{T}, U::Matrix{T}, H::Matrix{T}) where {T}
    common_iter!(HalleyUpdater(), X, U, H, alg.maxiter, alg.verbose, alg.tol)
end

struct HalleyUpdater <: PolarUpdater end


function update_U!(upd::HalleyUpdater, U::Matrix{T}) where {T}
    UtU = Array{T}(undef, size(U))
    mul!(UtU, transpose(U), U)
    copyto!(U, U * (3*I + UtU)* inv(I + 3*UtU))
end

#
# QR-based Dynamically Weighted Halley (QDWH) algorithm
#
# Reference: Optimizing Halley's iteration for computing the matrix 
#            polar decomposition, Yuji Nakatsukasa, Zhaojun Bai and 
#            Francois Gygi, SIAM, J. Mat. Anal. Vol. 31, Num 5 (2010)
#            pp. 2700-2720
#
# Limitations: 1. the QDWH method should support m > n matrix itself in the future.
#              2. the computing of alpha and L in solve! can be improved by using
#                 norm and condition number estimate, which are currently 
#                 not available in Julia.
#
mutable struct QDWHAlg{T} <: PolarAlg
    maxiter::Int
    verbose::Bool
    piv::Bool       # whether to pivot  
    tol::T  
    
    function QDWHAlg{T}( ;maxiter::Integer=100,
                     verbose::Bool=false,
                     piv::Bool=true, 
                     tol::Real=cbrt(eps(T))) where T
        maxiter > 1 || error("maxiter must be greater than 1.")
        tol > 0 || error("tol must be positive.")
        
        new{T}(maxiter,
            verbose,
            piv,
            tol)
    end
end


function solve!(alg::QDWHAlg,
                X::Matrix{T}, U::Matrix{T}, H::Matrix{T}) where {T}
    # alpha is an estimate of the largest singular value of the
    # original matrix
    X_temp = Array{T}(undef, size(X))
    copyto!(X_temp, X)

    n = size(X_temp, 1)
    alpha = norm(X_temp)   
    for i in length(X_temp)
        X_temp[i] /= alpha # form X0
    end

    # L is a lower bound for the smallest singular value of X0
    smin_est = opnorm(X_temp, 1)/cond(X_temp, 1)
    L  = smin_est/convert(T, sqrt(n))

    common_iter!(QDWHUpdater(alg.piv, L), X, U, H, alg.maxiter, alg.verbose, alg.tol)
end

mutable struct QDWHUpdater{T} <: PolarUpdater
    piv::Bool   # whether to pivot QR factorization
    L::T  # a lower bound for the smallest singluar value of each update matrix U
end


function update_U!(upd::QDWHUpdater, U::Matrix{T}) where {T}
    piv = upd.piv
    L = upd.L
    m, n = size(U)
    B = Array{T}(undef, m+n, n)
    Q1 = Array{T}(undef, n, n)
    Q2 = Array{T}(undef, n, n)
    # Compute paramters L, a, b, c
    L2 = L^2
    dd = try
        (4 * (1 - L2)/L2^2)^(1/3)
        catch 
        (complex(4 * (1 -L2)/L2^2, 0))^(1/3)
    end
    sqd = sqrt(1+dd)
    a = sqd + 0.5 * sqrt(8 - 4 * dd + 8 * (2 - L2)/(L2 * sqd))
    a = real(a)
    b = (a - 1)^2 / 4
    c = a + b - 1
    
    # update L
    upd.L = L * (a + b * L2)/(1 + c * L2)
    
    copyto!(B, [sqrt(c)*U; Matrix(one(T)*I,n,n)])
    if piv 
        F = qr(B, Val(true))
    else
        F = qr(B)
    end
    copyto!(Q1, Matrix(F.Q)[1:m, :])
    copyto!(Q2, Matrix(F.Q)[m+1:end, :])
    copyto!(U, b / c * U + (a - b / c) / sqrt(c) * Q1 * Q2')
    
end

