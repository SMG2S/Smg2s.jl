using Smg2s, SparseArrays, LinearAlgebra, Test

size = 8
offset = 3
vec=[1; 1; 0; 1; 1; 1; 0]
nilp = Nilp(vec, offset, size)
#nilp.nilpMat
M=sprandn(size, size, 0.75)

@test NilpxM(M, nilp) == nilp.nilpMat*M
@test MxNilp(M, nilp) == M*nilp.nilpMat
