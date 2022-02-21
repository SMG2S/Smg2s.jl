using SMG2S, Test, SparseArrays

@testset "Nilpotent Matrix" begin
    @testset "Nilp 1" begin
        @test SMG2S.Nilp(2,8).nilpMat * SMG2S.Nilp(2,8).nilpMat * SMG2S.Nilp(2,8).nilpMat == spzeros(8, 8)
    end
    @testset "Nilp 2" begin
        vec=[1; 1; 0; 1; 1; 1; 0]
        nilp = SMG2S.Nilp(vec, 8)
        @test nilp.nilpMat *nilp.nilpMat * nilp.nilpMat * nilp.nilpMat == spzeros(8, 8)
    end
    @testset "Nilp 3" begin
        nilpMatrix = sparse([0 1 1 0; 0 0 1 1 ; 0 0 0 1 ; 0 0 0 0])
        nilp=SMG2S.Nilp(nilpMatrix,4)
        @test nilp.nilpMat *nilp.nilpMat * nilp.nilpMat * nilp.nilpMat == spzeros(4, 4)
    end
    @testset "Nilp 4" begin
        vec=[1; 1; 0; 1; 1; 1; 0]
        nilp = SMG2S.Nilp(vec, 3, 8)
        @test nilp.nilpMat *nilp.nilpMat * nilp.nilpMat == spzeros(8, 8)
    end
    @testset "Nilp 5" begin
        nilp = SMG2S.Nilp(4, 3, 15)
        @test nilp.nilpMat * nilp.nilpMat * nilp.nilpMat * nilp.nilpMat * nilp.nilpMat == spzeros(15, 15)
    end
end
