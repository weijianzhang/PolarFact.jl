@testset "Newton" begin
    for n in [1,3,10]

A = rand(n,n)

r = polarfact(A, alg = :newton);

# Test unitary matrix U

U = r.U
H = r.H

@test U'*U ≈ Matrix(I,n,n) atol=1e-7

# Test Hermitian positive semifefinite matrix H

@test ishermitian(H)
 
for i in eigvals(H)
    @test i >= 0.
end

@test A ≈ U*H

    end
end
