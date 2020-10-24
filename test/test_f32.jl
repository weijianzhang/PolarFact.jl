@testset "Float32" begin
    for n in [1,3,10]

A = Array{Float32}(undef, n, n)
copyto!(A, rand(n,n))

r1 = polarfact(A, alg =:newton);

r2 = polarfact(A, alg =:qdwh);

r3 = polarfact(A, alg =:halley);

r4 = polarfact(A, alg =:svd);

r5 = polarfact(A, alg =:hybrid);

@test r1.U ≈ r2.U
@test r2.U ≈ r3.U
@test r3.U ≈ r4.U
@test r4.U ≈ r5.U

    end
end


