using Smg2s
using SparseArrays
using MatrixMarket
using LinearAlgebra

n = 10
diag_l = -5
diag_u = -3
nbone = 2

function f(idx::Integer, n = n)
    return 10 * cos(((idx-1) / 2) * 2 * π / n) + 5
end

spec = zeros(ComplexF64, n)
Spectrum!(spec, f, n)

genMat = nonsym(nbone, n, diag_l, diag_u, spec)

MatrixMarket.mmwrite("genMat.mtx", genMat)
