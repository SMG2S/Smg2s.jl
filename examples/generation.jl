using Smg2s
using Test
using SparseArrays
using Plots

using LinearAlgebra
using Random

n = 1000
diag_l = -800
diag_u = -790
nbOne = 10

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

nilpMatrix = sparse(diagm(100=>nilpvec[1:n-100]))

#@info nilpvec
#nilp = Nilp(nilpMatrix, n)

nilp = Nilp(nbOne, 100, n)

## No hermitian case
function f1(idx::Integer, n = n)
    return cos((idx-1) * 2 * π / n) + 1 * rand() + 10 + 10 * sin((idx-1) * 2 * π / n) * im
end

spec1 = zeros(ComplexF64, n)
Spectrum!(spec1, f1, n)

Am1 = spzeros(ComplexF64, n, n)
initMat!(Am1, diag_l, diag_u, n; scale=0.1, sparsity=0.5)

@time genMat1 = nonherm(n, diag_l, diag_u, spec1, Am1, nilp)
@info "sparsity = " nnz(genMat1) / (n * n)
genspec1 = eigen!(Matrix(genMat1)).values

p1=scatter(real(spec1), imag(spec1), markercolor = :green, markern = 7, label="given spectrum")
p1=scatter!(real(genspec1), imag(genspec1), markercolor = :red, markern = 5, label="spectrum of generated matrix")

## Non symmetric 1
function f2(idx::Integer, n = n)
    rnd = rand()
    if idx % 2 == 0
        return 10 * cos((idx) * 2 * π / n) + 5 * rnd + 5 * sin((idx) * 2 * π / n) * im
    else
        return 10 * cos((idx+1) * 2 * π / n) + 5 * rnd - 5 * sin((idx+1) * 2 * π / n) * im
    end
end

spec2 = zeros(ComplexF64, n)
Spectrum!(spec2, f2, n)


Am2 = spzeros(Float64, n, n)
initMat!(Am2, diag_l, diag_u, n; scale=0.1, sparsity=0.5)
@time genMat2 = nonsym(n, diag_l, diag_u, spec2, Am2, nilp)

@info "sparsity = " nnz(genMat2) / (n * n)

genspec2 = eigen!(Matrix(genMat2)).values
p2=scatter(real(spec2), imag(spec2), markercolor = :green, markern = 7, label="given spectrum")
p2=scatter!(real(genspec2), imag(genspec2), markercolor = :red, markern = 5, label="spectrum of generated matrix")


## Non symmetric 2
function f3(idx::Integer, n = n)
    return 10 * cos(((idx-1) / 2) * 2 * π / n) + 5
end

spec3 = zeros(ComplexF64, n)
Spectrum!(spec3, f3, n)


Am3 = spzeros(Float64, n, n)
initMat!(Am3, diag_l, diag_u, n; scale=0.1, sparsity=.5)
@time genMat3 = nonsym(n, diag_l, diag_u, spec3, Am3, nilp)

@info "sparsity = " nnz(genMat3) / (n * n)

genspec3 = eigen!(Matrix(genMat3)).values
p3=scatter(real(spec3), imag(spec3), ylims=(-1,1), markercolor = :green, markern = 7, label="given spectrum")
p3=scatter!(real(genspec3), imag(genspec3), ylims=(-1,1), markercolor = :red, markern = 5, label="spectrum of generated matrix")

p4=spy(dropzeros(genMat3), legend = :none)

l = @layout [a  b;c d]
plot(p1, p2, p3, p4, layout = l)

#savefig("fig/example")
