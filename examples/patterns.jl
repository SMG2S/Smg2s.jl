using SMG2S
using Test
using SparseArrays
using Plots

using LinearAlgebra
using Random

size = 10000
diag_l = -9900
diag_u = -800
nbOne = 70
nilpdiag = 3500
sp=0.0002
figname = "fig/pattern12"

nilpvec = zeros(Int64, size-1)

cnt = 0

while cnt < size-1
    l₁ = rand(1:nbOne)
    l₀ = rand(1:nbOne)

    for i = cnt+1:l₁+cnt
        if i <= size-1
            nilpvec[i] = 1
        end
    end

    global cnt += l₀ + l₁

end

nilpMatrix = sparse(diagm(nilpdiag=>nilpvec[1:size-nilpdiag]))

nilp = Nilp(nilpMatrix, size)

function f(idx::Integer, size = size)
    return 10 * cos(((idx-1) / 2) * 2 * π / size) + 5
end

spec = zeros(ComplexF64, size)
Spectrum!(spec, f, size)

Am = spzeros(Float64, size, size)
initMat!(Am, diag_l, diag_u, size; scale=0.1, sparsity=sp)
genMat = nonsym(size, diag_l, diag_u, spec, Am, nilp)
sparsity = nnz(genMat) / (size * size)
@info "sparsity=" sparsity
spy(genMat, legend = :none, title="size:$size, diag_l:$diag_l, diag_u:$diag_u, nbOne:$nbOne, nilpdiag: $nilpdiag, sp: $sp",titlefont = font(8))

savefig(figname)

#savefig("fig/example")
