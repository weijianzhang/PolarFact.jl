using PolarFact
using Test, LinearAlgebra
using Random

Random.seed!(1103)

include("test_newton.jl")
include("test_halley.jl")
include("test_svd.jl")
include("test_hybrid.jl")
include("test_f32.jl") # test Float32 case
