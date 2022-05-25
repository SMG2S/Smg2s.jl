using Smg2s, SparseArrays, LinearAlgebra, Test

dim = 5000
nzeros = 10
nilp = Nilp(10, dim)
#nilp.nilpMat
M=sprandn(dim, dim, 0.25)

@test MNilpM(M, nilp) == M*NilpMat(nilp) - NilpMat(nilp)*M

@time MNilpM(M, nilp)
@time M*NilpMat(nilp) - NilpMat(nilp)*M

@info "done"
