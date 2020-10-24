@testset "Halley" begin
    for n in [1,3,10]

A = rand(n, n)

r = polarfact(A, alg = :halley);


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

##########################################################
@testset "QDWH" begin
    for n in [1,3,10]

m = rand(1:10)

B = rand(m, m)

r2 = polarfact(B, alg = :qdwh);

# Test unitary matrix U

U2 = r2.U
H2 = r2.H

@test U2'*U2 ≈ Matrix(I,m,m) atol=1e-7

# Test Hermitian positive semifefinite matrix H

@test issymmetric(H2)
 
for i in eigvals(H2)
    @test i >= 0.
end

@test B ≈ U2*H2


    end
end
