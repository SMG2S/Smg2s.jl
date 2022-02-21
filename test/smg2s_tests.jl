using Smg2s, Test, SparseArrays, LinearAlgebra

#Define size of matrix to be generated
n = 20

#maximum number of continous `1` in the diagonal of generated nilpotent matrix
nbOne = 3
#offsets of lower diagonals between which the matrix will be initially filled
diag_l = -10
diag_u = -5


@testset "non-Herm SMG2S" begin
    @testset "Generation 1" begin
        #a function to generate the user-provided spectrum, `idx` is the indexing of eigenvalues
        function f(idx::Integer, n = n)
            return 10 * idx + 1 + idx * im
        end
        spec = zeros(ComplexF64, n)
        Spectrum!(spec, f, n)
        A = nonherm(nbOne, n, diag_l, diag_u, spec)
        @test spec ≈ eigen!(Matrix(A)).values atol=0.1
    end
end

@testset "non-Sym SMG2S" begin
    @testset "Generation 1" begin
        #a function to generate the user-provided spectrum, `idx` is the indexing of eigenvalues
        function f(idx::Integer, n = n)
            if idx % 2  == 1
                return  ((idx + 1) / 2)  + 1 - 0.1 * ((idx + 1) / 2) * im
            else
                return  (idx / 2) + 1 + 0.1 * (idx / 2) * im
            end
        end
        spec = zeros(ComplexF64, n)
        Spectrum!(spec, f, n)
        A = nonsym(nbOne, n, diag_l, diag_u, spec)
        @test spec ≈ eigen!(Matrix(A)).values atol=0.5
    end
end
