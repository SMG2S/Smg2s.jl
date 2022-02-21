using Smg2s
using Test
using SparseArrays
using Plots

using LinearAlgebra
using Random

n = 10000
diag_l = -9900
diag_u = -800
nbOne = 70
nilpdiag = 3500
sp=0.0002
figname = "fig/pattern12"

nilpvec = zeros(Int64, n-1)

cnt = 0

while cnt < n-1
    l₁ = rand(1:nbOne)
    l₀ = rand(1:nbOne)

    for i = cnt+1:l₁+cnt
        if i <= n-1
            nilpvec[i] = 1
        end
    end

    global cnt += l₀ + l₁

end

nilpMatrix = sparse(diagm(nilpdiag=>nilpvec[1:n-nilpdiag]))

nilp = Nilp(nilpMatrix, n)

function f(idx::Integer, n = n)
    return 10 * cos(((idx-1) / 2) * 2 * π / n) + 5
end

spec = zeros(ComplexF64, n)
Spectrum!(spec, f, n)

Am = spzeros(Float64, n, n)
initMat!(Am, diag_l, diag_u, n; scale=0.1, sparsity=sp)
genMat = nonsym(n, diag_l, diag_u, spec, Am, nilp)
sparsity = nnz(genMat) / (n * n)
@info "sparsity=" sparsity
spy(genMat, legend = :none, title="size:$n, diag_l:$diag_l, diag_u:$diag_u, nbOne:$nbOne, nilpdiag: $nilpdiag, sp: $sp",titlefont = font(8))
