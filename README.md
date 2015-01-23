# PolarFact

[![Build Status](https://travis-ci.org/weijianzhang/PolarFact.jl.svg?branch=master)](https://travis-ci.org/weijianzhang/PolarFact.jl)
| Julia 0.3 [![PolarFact](http://pkg.julialang.org/badges/PolarFact_release.svg)](http://pkg.julialang.org/?pkg=PolarFact&ver=release)
| Julia 0.4 [![PolarFact](http://pkg.julialang.org/badges/PolarFact_nightly.svg)](http://pkg.julialang.org/?pkg=PolarFact&ver=nightly)

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
nearest orthonormal matrix to ``A`` in the Frobenius norm [1] (Sec. 8.1). 

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
backward stable [6]. The QDWH is a cubic-rate convergent method.
It is backward stable under the assumption that column pivoting and
either row pivoting or row sorting are used in the QR factorization [6]. 
Without scaling, both type of methods can be slow when the matrix is
ill-conditioned.

On many modern computers, matrix multiplication can be performed
very efficiently. The Newton Schulz method requires two matrix
multiplication while the (scaled) Newton method requires one matrix
inversion. Thus the hybrid Newton is more efficient if matrix
multiplication is 2 times faster than the matrix inversion [4].

Comparing to the SVD approach, the iterative algorithms are much more
efficient when the matrix is nearly unitary (arises in aerospace
applications). 

[6] Yuji Nakatsukasa and Nicholas J. Higham, Backward stability of
iterations for computing the polar decomposition, SIAM, J.
Matrix Anal. Appl. Vol. 33, No. 2, pp. 460-479. 


## Interface

The package provides a high-level function ``polarfact``:

```julia
	polarfact(A; alg, maxiter, tol, verbose)
```

The meaning of the arguments:

- ``A`` : the input matrix of type ``Matrix{Float64}``.

- ``alg``: a symbol that indicates the factorization algorithm (default = ``:newton``).

	This argument accepts the following values:

	- ``:newton``: scaled Newton's method
	- ``:qdwh``: the QR-based Dynamically weighted Halley (QDWH) method
	- ``:halley``: Halley's method
	- ``:schulz``: the Newton Schulz method
	- ``:hybrid``: a hybrid Newton method 
	- ``:svd``: the SVD method

- ``maxiter``: maximum number of iterations (default = ``100``).

- ``tol`` :  tolerance (default = ``1.0e-6``).

- ``verbose`` : whether to show procedural information (default = ``false``).

*Note:* ``maxiter``, ``tol`` and ``verbose`` are not used for the
SVD method.

The output has type ``PolarFact.Result``, which is defined as 

```
	immutable Result
		U::Matrix{Float64}               # unitary factor
		H::Matrix{Float64}               # Hermitian positive semidefinite factor
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

julia> r = polarfact(A, verbose=true);
Iter.    Rel. err.        Obj.         
    1     8.031554e-01     9.218344e+00
    2     4.931430e-01     4.609596e-01
    3     1.073877e-01     8.672101e-03
    4     2.239146e-03     8.894508e-06
    5     2.256007e-06     1.782681e-11
    6     4.474254e-12     6.522560e-16

julia> r.U
6x6 Array{Float64,2}:
  0.51827    0.315836    0.563476      0.396524   -0.266376   -0.293174 
  0.674541  -0.164589   -0.000399023  -0.266933   -0.0162729   0.668122 
 -0.401425  -0.32354     0.654859      0.287182    0.158807    0.444575 
  0.206237  -0.0348974  -0.392979      0.72679     0.515961    0.0858893
  0.103665   0.291805    0.303866     -0.393481    0.792383   -0.170502 
 -0.248935   0.825901   -0.082991      0.0774023  -0.0973267   0.483289 

julia> r.niters
6

julia> using MatrixDepot

julia> B = matrixdepot("hilb", 6) # try Hilbert matrix
6x6 Array{Float64,2}:
 1.0       0.5       0.333333  0.25      0.2       0.166667 
 0.5       0.333333  0.25      0.2       0.166667  0.142857 
 0.333333  0.25      0.2       0.166667  0.142857  0.125    
 0.25      0.2       0.166667  0.142857  0.125     0.111111 
 0.2       0.166667  0.142857  0.125     0.111111  0.1      
 0.166667  0.142857  0.125     0.111111  0.1       0.0909091

julia> r = polarfact(B, alg = :halley, verbose=true);
Iter.    Rel. err.        Obj.         
    1     4.980074e-01     1.648164e+00
    2     4.202594e-01     1.676818e+00
    3     3.111776e-01     1.642458e+00
    4     5.572934e-01     1.607827e+00
    5     1.490292e-01     1.585184e+00
    6     2.150641e-01     1.518917e+00
    7     3.984436e-01     1.466951e+00
    8     1.123055e-01     1.452383e+00
    9     1.494373e-01     1.423552e+00
   10     3.788419e-01     1.256390e+00
   11     3.466910e-01     1.278821e+00
   12     3.509951e-02     1.281831e+00
   13     1.071601e-01     1.248502e+00
   14     2.977220e-01     9.945915e-01
   15     4.693540e-01     2.107683e-01
   16     1.063725e-01     4.599108e-04
   17     2.299602e-04     3.678240e-12
   18     1.839120e-12     2.220446e-16

julia> r.niters
18

julia> r = polarfact(B, alg = :newton, verbose=true);
Iter.    Rel. err.        Obj.         
    1     1.236974e+03     6.297632e+06
    2     9.907912e-01     8.240846e+02
    3     9.407495e-01     4.673527e+00
    4     5.931747e-01     7.838535e-02
    5     3.697098e-02     3.702364e-04
    6     1.850680e-04     7.358231e-09
    7     3.679115e-09     7.669935e-18

julia> r.niters
7
```
## Acknowledgements

The design of the package is inspired by [NMF.jl](https://github.com/JuliaStats/NMF.jl).


