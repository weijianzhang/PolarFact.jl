# PolarFact

[![Build Status](https://travis-ci.com/weijianzhang/PolarFact.jl.svg?branch=master)](https://travis-ci.org/weijianzhang/PolarFact.jl)

A Julia package for the matrix polar decomposition.

## Install

To install the release version, type

```julia
julia> Pkg.add("PolarFact")
```

## Overview 

Every ``m-by-n`` matrix ``A`` has a polar decomposition ``A=UH``,
where the ``m-by-n`` matrix ``U`` has orthonormal columns if ``m>n``
or orthonormal rows if ``m<n`` and the ``n-by-n`` matrix ``H`` is
Hermitian positive semidefinite. For a square matrix ``A``, ``H`` is
unique. If in addition, ``A`` is nonsingular, then ``H`` is positive
definite and ``U`` is unique.

The polar decomposition is closely related to the singular value
decomposition (SVD). In particular, if ``A = P*S*Q'`` is a singular
value decomposition of A, then ``U = P*Q'`` and ``H = Q*S*Q'`` are the
corresponding polar factors. The orthonormal polar factor ``U`` is the
nearest orthonormal matrix to ``A`` in the Frobenius norm [1, Sec. 8.1]. 

[1] Nicholas J. Higham, Functions of Matrices: Theory and Computation,
SIAM, Philadelphia, PA, USA, 2008.

This package provides the following algorithms for computing matrix
polar decomposition:

* (Scaled) Newton's method

	Reference:
	[2] Nicholas J. Higham, Computing the Polar Decomposition ---with Applications,
	SIAM J. Sci. Statist. Comput. Vol. 7, Num 4 (1986) pp. 1160-1174.
	
* the Newton Schulz method 
  
    This method can only apply to matrix ``A`` such that ``norm(A) < sqrt(3)``.

	Reference:
	[3] Günther Schulz, Iterative Berechnung der reziproken Matrix, Z. Angew.
	Math. Mech.,13:57-59, (1933) pp. 114, 181.

* a hybrid Newton method

	Start with (scaled) Newton's method and switch to Newton-Schulz method
	when convergence is guaranteed.

	Reference:
	[4] Nicholas J. Higham and Robert S. Schreiber, Fast Polar
	Decomposition of an arbitrary matrix, SIAM, J. Sci. Statist. Comput.
	Vol. 11, No. 4 (1990) pp. 648-655

* Halley's method

	Reference:
	[5] Y. Nakatsukasa, Z. Bai and F. Gygi, Optimizing Halley's iteration 
	for computing the matrix polar decomposition, SIAM, J. Mat. Anal. 
	Vol. 31, Num 5 (2010) pp. 2700-2720. 

* the QR-based Dynamically weighted Halley (QDWH) method [5]  

* the SVD method

### Comments on Usage

The scaled Newton iteration is a well known and effective method for
computing the polar decomposition. It converges quadratically and is
backward stable under the assumption that the matrix inverses are
computed in a mixed backward forward stable way [6]. The QDWH is a
cubic-rate convergent method.  It is backward stable under the
assumption that column pivoting and either row pivoting or row sorting
are used in the QR factorization [6].  Without scaling, both type of
methods can be slow when the matrix is ill-conditioned.

On many modern computers, matrix multiplication can be performed
very efficiently. The Newton Schulz method requires two matrix
multiplication while the (scaled) Newton method requires one matrix
inversion. Thus the hybrid Newton is more efficient if matrix
multiplication is 2 times faster than the matrix inversion [4].

Comparing to the SVD approach, the iterative algorithms are much more
efficient when the matrix is nearly unitary (as in applications, for
example, where a time-dependent matrix drifts from orthogonality due
to rounding errors or other errors).

[6] Yuji Nakatsukasa and Nicholas J. Higham, Backward stability of
iterations for computing the polar decomposition, SIAM, J.
Matrix Anal. Appl. Vol. 33, No. 2, pp. 460-479. 


## Interface

The package provides a high-level function ``polarfact``:

```julia
	polarfact(A; alg, maxiter, tol, verbose)
```

The meaning of the arguments:

- ``A`` : the input matrix of type ``Matrix{T}``, where ``T`` could be
          ``Float16``, ``Float32``, ``Float64``, ``BigFloat``. 

- ``alg``: a symbol that indicates the factorization algorithm (default = ``:newton``).

	This argument accepts the following values:

	- ``:newton``: scaled Newton's method
	- ``:qdwh``: the QR-based Dynamically weighted Halley (QDWH) method
	- ``:halley``: Halley's method
	- ``:schulz``: the Newton Schulz method
	- ``:hybrid``: a hybrid Newton method 
	- ``:svd``: the SVD method

- ``maxiter``: maximum number of iterations (default = ``100``).

- ``tol`` :  tolerance (default = ``cbrt(eps(T))``).

- ``verbose`` : whether to show procedural information (default = ``false``), where
               ``Iter`` is the number of iterations, ``Rel. err.`` is equal to
			   ``norm(preU - U, Inf) / norm(preU, Inf)`` and ``Obj.`` is equal to
			   ``norm(I - U'*U, Inf)``. 

*Note:* ``maxiter``, ``tol`` and ``verbose`` are not used for the
SVD method.

The output has type ``PolarFact.Result``, which is defined as 

```
	immutable Result{T}
		U::Matrix{T}               # unitary factor
		H::Matrix{T}               # Hermitian positive semidefinite factor
		niters::Union(Int, Nothing)      # number of iterations or Nothing
		converged::Union(Bool, Nothing)  # whether the algorithm converges or Nothing
	end
```

*Note:* ``niters`` and ``converged`` are of type ``Nothing`` for the
SVD method. 

## Examples

```julia
julia> using PolarFact

julia> A = rand(6,6);

julia> r = polarfact(A);

julia> r.U
6x6 Array{Float64,2}:
  0.78067    -0.0694445   0.470076    -0.292781   -0.0423338   0.277934 
 -0.0984991   0.200495   -0.350106    -0.167351    0.394071    0.802638 
  0.235931    0.376468    0.00966701   0.0259734   0.78062    -0.438716 
  0.38059    -0.347928   -0.256805     0.79899     0.105905    0.136195 
  0.21592     0.816347   -0.136636     0.261595   -0.445576    0.0362852
 -0.365675    0.160542    0.756137     0.422827    0.154252    0.257271 

julia> r.niters
6

julia> using MatrixDepot  # a test matrix collection

julia> A = matrixdepot("randsvd", 20, 10^15);  # test a very ill conditioned random matrix 

julia> r = polarfact(A, alg = :newton, verbose = true);
Iter.    Rel. err.        Obj.         
    1     1.548278e+07     4.804543e+14
    2     9.999002e-01     5.902800e+06
    3     9.925194e-01     7.247316e+02
    4     9.317962e-01     9.261149e+00
    5     7.226745e-01     3.409441e-01
    6     1.392861e-01     2.550612e-03
    7     1.272341e-03     2.760504e-07
    8     1.378506e-07     1.062337e-14

julia> r = polarfact(A, alg = :qdwh, verbose = true);
Iter.    Rel. err.        Obj.         
    1     1.018823e+00     2.294023e+00
    2     9.492166e-01     2.113416e+00
    3     6.766440e-01     7.363896e-01
    4     1.337401e-01     8.038009e-04
    5     1.208881e-04     4.126278e-13
    6     6.184049e-14     2.242130e-15

julia> r = polarfact(A, alg = :halley);

julia> r.niters
34

julia> A = matrixdepot("hadamard", Float32, 8) # test Float32 Hadamard matrix
8x8 Array{Float32,2}:
 1.0   1.0   1.0   1.0   1.0   1.0   1.0   1.0
 1.0  -1.0   1.0  -1.0   1.0  -1.0   1.0  -1.0
 1.0   1.0  -1.0  -1.0   1.0   1.0  -1.0  -1.0
 1.0  -1.0  -1.0   1.0   1.0  -1.0  -1.0   1.0
 1.0   1.0   1.0   1.0  -1.0  -1.0  -1.0  -1.0
 1.0  -1.0   1.0  -1.0  -1.0   1.0  -1.0   1.0
 1.0   1.0  -1.0  -1.0  -1.0  -1.0   1.0   1.0
 1.0  -1.0  -1.0   1.0  -1.0   1.0   1.0  -1.0

julia> r = polarfact(A)
Result{Float32}(8x8 Array{Float32,2}:
 0.353553   0.353553   0.353553   0.353553  …   0.353553   0.353553   0.353553
 0.353553  -0.353553   0.353553  -0.353553     -0.353553   0.353553  -0.353553
 0.353553   0.353553  -0.353553  -0.353553      0.353553  -0.353553  -0.353553
 0.353553  -0.353553  -0.353553   0.353553     -0.353553  -0.353553   0.353553
 0.353553   0.353553   0.353553   0.353553     -0.353553  -0.353553  -0.353553
 0.353553  -0.353553   0.353553  -0.353553  …   0.353553  -0.353553   0.353553
 0.353553   0.353553  -0.353553  -0.353553     -0.353553   0.353553   0.353553
 0.353553  -0.353553  -0.353553   0.353553      0.353553   0.353553  -0.353553,8x8 Array{Float32,2}:
  2.82843     -1.49012e-8  -5.96046e-8  …  -1.49012e-8  -1.19209e-7  -1.49012e-8
 -1.49012e-8   2.82843      5.96046e-8      5.96046e-8  -4.47035e-8  -7.45058e-8
 -5.96046e-8   5.96046e-8   2.82843         0.0          8.9407e-8   -8.9407e-8 
 -7.45058e-8   1.49012e-8   2.98023e-8      1.3411e-7    1.3411e-7    1.49012e-7
 -2.98023e-8   4.47035e-8   0.0             1.49012e-8   2.98023e-8  -7.45058e-8
 -1.49012e-8   5.96046e-8   0.0         …   2.82843      1.63913e-7   7.45058e-8
 -1.19209e-7  -4.47035e-8   8.9407e-8       1.63913e-7   2.82843      7.45058e-8
 -1.49012e-8  -7.45058e-8  -8.9407e-8       7.45058e-8   7.45058e-8   2.82843   ,2,true)

```

## Acknowledgements

The design of the package is inspired by [NMF.jl](https://github.com/JuliaStats/NMF.jl).


