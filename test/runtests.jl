using SMG2S
using Test
using SparseArrays
using Plots

using LinearAlgebra
using Random

size = 1000
diag_l = -10
diag_u = -2
nbOne = 8

function f(idx::Integer, size = size)
    return cos((idx-1) * 2 * π / size) + 1 * rand() + 10 + 10 * sin((idx-1) * 2 * π / size) * im
end

spec = zeros(ComplexF64, size)
Spectrum!(spec, f, size)

Am = spzeros(ComplexF64, size, size)
initMat!(Am, diag_l, diag_u, size; scale=0.1, sparsity=0.5)

@time genMat = nonherm(nbOne, size, diag_l, diag_u, spec, Am)
@info "sparsity = " nnz(genMat) / (size * size)
genspec = eigen!(Matrix(genMat)).values

scatter(real(spec), imag(spec), markercolor = :green, markersize = 7, label="given spectrum")
scatter!(real(genspec), imag(genspec), markercolor = :red, label="spectrum of generated matrix")
