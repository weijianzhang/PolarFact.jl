@testset "SVD" begin
    for n in [1,3,10]

m = n + rand(0:5)
A = rand(m, n)

r = polarfact(A, alg = :svd);


# Test unitary matrix U

U = r.U
H = r.H

@test U'*U ≈ Matrix(I,n,n) atol=1e-7

# Test Hermitian positive semifefinite matrix H

@test issymmetric(H)
 
for i in eigvals(H)
    @test i >= 0.
end

@test A ≈ U*H

    end
end

