# common types or functions

# the result type

struct Result{T}
    U::Matrix{T}
    H::Matrix{T}
    niters::Integer
    converged::Bool
    
    function Result(U::Matrix{T}, H::Matrix{T}, 
                    niter::Integer, 
                    converged::Bool) where {T}
        size(U, 2) == size(H, 1) || 
               throw(DimensionMismatch("Inner dimension of U and H mismatch."))
        new{T}(U, H, niter, converged)
    end
end

struct SVDResult{T}
    U::Matrix{T}
    H::Matrix{T}
    function SVDResult(U::Matrix{T}, H::Matrix{T}) where {T}
        size(U, 2) == size(H, 1) || 
               throw(DimensionMismatch("Inner dimension of U and H mismatch."))
        new{T}(U, H)
    end
end

# niters and converged are nothing for the SVD method

# the objective type

struct Objective{T}
    absolute::T  # absolute error
    relative::T  # relative error

    function Objective(absolute::T, relative::T) where {T <: Real}
        absolute >= 0 || error("absolute must be non-negative.")
        relative >= 0 || error("relative must be non-negative.")

        new{T}(absolute, relative)
    end
end


function evaluate_objv(preU::Matrix{T}, U::Matrix{T}) where {T}
    rel = opnorm(preU - U, Inf) / opnorm(preU, Inf)
    abs = opnorm(I - U'*U, Inf)
    return Objective{T}(rel, abs)
end


abstract type PolarUpdater end

abstract type PolarAlg end

# common algorithm skeleton for iterative updating methods

function common_iter!(updater::PolarUpdater,
                         X::Matrix{T}, 
                         U::Matrix{T}, 
                         H::Matrix{T},
                         maxiter::Int,
                         verbose::Bool,
                         tol::T) where {T}
    
    preU = Array{T}(undef, size(X))
    copyto!(U, X)
    converged = false
    t = 0
    if verbose
        @printf("%-5s    %-13s    %-13s\n", "Iter.", "Rel. err.",  "Obj.")
    end

    while !converged && t < maxiter
        t += 1
        copyto!(preU, U)
        update_U!(updater, U)
               
        # determine convergence
        diff = norm(preU - U)
        if diff < tol
            converged = true
        end
        
        # display infomation
        if verbose
            objv = evaluate_objv(preU, U)
            @printf("%5d    %13.6e    %13.6e\n", 
                    t, objv.absolute, objv.relative)
        end
    end

    # compute H
    mul!(H, U', X)
    H = (1/T(2)) * (H + H')
    return Result(U, H, t, converged)
end

# Scaling iterative algorithm
function common_iter_scal!(updater::PolarUpdater,
                              X::Matrix{T}, 
                              U::Matrix{T}, 
                              H::Matrix{T},
                              maxiter::Int,
                              verbose::Bool,
                              tol::T) where {T}
    
    preU = Array{T}(undef, size(X))
    copyto!(U, X)
    converged = false
    t = 0
    if verbose
        @printf("%-5s    %-13s    %-13s\n", "Iter.", "Rel. err.",  "Obj.")
    end

    while !converged && t < maxiter
        t += 1
        copyto!(preU, U)
        update_U!(updater, U)
        
        # determine convergence
        diff = norm(preU - U)
        if diff < tol
            converged = true
        end
        
        # determine scaling
        reldiff = diff/norm(U) # relative error
        if updater.scale && (reldiff < updater.scale_tol)
            updater.scale = false
        end

        # display infomation
        if verbose
            objv = evaluate_objv(preU, U)
            @printf("%5d    %13.6e    %13.6e\n", 
                    t, objv.absolute, objv.relative)
        end
    end

    # compute H
    mul!(H, U', X)
    H = (1/T(2)) * (H + H')
    return Result(U, H, t, converged)
end



# Hybrid iteration algorithm 
function common_iter_hybr!(updater1::PolarUpdater,
                              updater2::PolarUpdater,
                              X::Matrix{T}, 
                              U::Matrix{T}, 
                              H::Matrix{T},
                              maxiter::Int,
                              verbose::Bool,
                              tol::T,
                              theta::T) where {T} # theta is the switch parameter

    preU = Array{T}(undef, size(X))
    copyto!(U, X)
    converged = false
    switched = false
    t = 0
    if verbose
        @printf("%-5s    %-13s    %-13s\n", "Iter.", "Rel. err.",  "Obj.")
    end

    while !converged && t < maxiter
        t += 1
        copyto!(preU, U)

        if switched
            update_U!(updater2, U)        
        else
            obj = opnorm(I - U'*U, 1)
            if obj > theta # theta is the switch parameter
                update_U!(updater1, U)
            else
                switched = true
                update_U!(updater2, U)
            end
        end

        # determine convergence
        diff = norm(preU - U)
        if diff < tol
            converged = true
        end
        
        # display infomation
        if verbose
            objv = evaluate_objv(preU, U)
            @printf("%5d    %13.6e    %13.6e\n", 
                    t, objv.absolute, objv.relative)
        end
    end

    # compute H
    mul!(H, U', X)
    H = (1/T(2)) * (H + H')
    return Result(U, H, t, converged)
end
