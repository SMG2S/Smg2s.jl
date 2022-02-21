using Smg2s, Test, SparseArrays

@testset "Nilpotent Matrix" begin
    @testset "Nilp 1" begin
        @test Nilp(2,8).nilpMat * Nilp(2,8).nilpMat * Nilp(2,8).nilpMat == spzeros(8, 8)
    end
    @testset "Nilp 2" begin
        vec=[1; 1; 0; 1; 1; 1; 0]
        nilp = Nilp(vec, 8)
        @test nilp.nilpMat *nilp.nilpMat * nilp.nilpMat * nilp.nilpMat == spzeros(8, 8)
    end
    @testset "Nilp 3" begin
        nilpMatrix = sparse([0 1 1 0; 0 0 1 1 ; 0 0 0 1 ; 0 0 0 0])
        nilp=Nilp(nilpMatrix,4)
        @test nilp.nilpMat *nilp.nilpMat * nilp.nilpMat * nilp.nilpMat == spzeros(4, 4)
    end
    @testset "Nilp 4" begin
        vec=[1; 1; 0; 1; 1; 1; 0]
        nilp = Nilp(vec, 3, 8)
        @test nilp.nilpMat *nilp.nilpMat * nilp.nilpMat == spzeros(8, 8)
    end
    @testset "Nilp 5" begin
        nilp = Nilp(4, 3, 15)
        @test nilp.nilpMat * nilp.nilpMat * nilp.nilpMat * nilp.nilpMat * nilp.nilpMat == spzeros(15, 15)
    end
end
