#
# Compute Polar Decomposition via SVD
#

struct SVDAlg <: PolarAlg end

function solve!(alg::SVDAlg, 
                X::Matrix{T}, U::Matrix{T}, H::Matrix{T}) where {T}

    F = svd(X, full = false)
    PQt = Array{T}(undef, size(U))
    mul!(PQt, F.U, transpose(F.V))
    copyto!(U, PQt)
    copyto!(H, F.V * Diagonal(F.S) * F.Vt)
    H = (1/T(2)) * (H + H')
    
    return SVDResult(U, H)
end
