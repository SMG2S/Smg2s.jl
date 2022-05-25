using Smg2s, Test, SparseArrays

@testset "Nilpotent Matrix" begin
    @testset "Nilp 1" begin
        @test NilpMat(Nilp(2,8)) * NilpMat(Nilp(2,8)) * NilpMat(Nilp(2,8)) == spzeros(8, 8)
    end
    @testset "Nilp 2" begin
        vec=[1; 1; 0; 1; 1; 1; 0]
        nilp = Nilp(vec, 8)
        @test NilpMat(nilp) *NilpMat(nilp) * NilpMat(nilp) * NilpMat(nilp) == spzeros(8, 8)
    end
    @testset "Nilp 3" begin
        vec=[1; 1; 0; 1; 1; 1; 0]
        nilp = Nilp(vec, 3, 8)
        @test NilpMat(nilp)*NilpMat(nilp) * NilpMat(nilp) == spzeros(8, 8)
    end
    @testset "Nilp 4" begin
        nilp = Nilp(4, 3, 15)
        @test NilpMat(nilp) * NilpMat(nilp) * NilpMat(nilp) * NilpMat(nilp) * NilpMat(nilp) == spzeros(15, 15)
    end
end
