using SMG2S
using SparseArrays
using MatrixMarket

function f(idx::Integer, size = size)
    return 10 * cos(((idx-1) / 2) * 2 * π / size) + 5
end

spec = zeros(ComplexF64, size)
Spectrum!(spec, f, size)

genMat = nonsym(size, diag_l, diag_u, spec)

MatrixMarket.mmwrite("genMat.mtx", genMat)
