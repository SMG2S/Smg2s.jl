using Smg2s, Test, SparseArrays, LinearAlgebra

@testset "Initialization of Matrix" begin
    @testset "Init 1" begin
        Am=spzeros(ComplexF64, 50, 50)
        initMat!(Am, -20, -10, 50, sparsity=0.1)
        @test nnz(Am) ≈ 38 atol=5
        @test UpperTriangular(Am)-diagm(diag(Am)) == zeros(50,50)
    end
    @testset "Init 2" begin
        Am = spzeros(Float64, 50, 50);
        initMat!(Am, -20, -10, 50, sparsity=0.05)
        @test nnz(Am) ≈ 19 atol=5
        @test UpperTriangular(Am)-diagm(diag(Am)) == zeros(50,50)
        @test diag(Am,-1) == zeros(49)
    end
end
