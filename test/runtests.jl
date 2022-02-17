using SMG2S
using Test
using SparseArrays
using Plots

using LinearAlgebra
using Random

size = 2000
diag_l = -120
diag_u = -5
nbOne = 40

nilpvec = zeros(Int64, size)

cnt = 0

while cnt <= size
    l₁ = rand(1:nbOne)
    l₀ = rand(1:nbOne)

    for i = cnt+1:l₁+cnt
        if i <= size
            nilpvec[i] = 1
        end
    end

    for i = cnt+l₁+1:l₁+l₀+cnt
        if i <= size
            nilpvec[i] = 0
        end
    end

    global cnt += l₀ + l₁

end

#@info nilpvec
nilp = Nilp(nilpvec, size)

## No hermitian case
function f1(idx::Integer, size = size)
    return cos((idx-1) * 2 * π / size) + 1 * rand() + 10 + 10 * sin((idx-1) * 2 * π / size) * im
end

spec1 = zeros(ComplexF64, size)
Spectrum!(spec1, f1, size)

Am1 = spzeros(ComplexF64, size, size)
initMat!(Am1, diag_l, diag_u, size; scale=0.1, sparsity=0.005)

@time genMat1 = nonherm(size, diag_l, diag_u, spec1, Am1, nilp)
@info "sparsity = " nnz(genMat1) / (size * size)
genspec1 = eigen!(Matrix(genMat1)).values

p1=scatter(real(spec1), imag(spec1), markercolor = :green, markersize = 7, label="given spectrum")
p1=scatter!(real(genspec1), imag(genspec1), markercolor = :red, markersize = 5, label="spectrum of generated matrix")

## Non symmetric 1
function f2(idx::Integer, size = size)
    if idx % 2 == 0
        return 10 * cos(((idx-1) / 2) * 2 * π / size) + 5 + 5 * sin(((idx-1) / 2) * 2 * π / size) * im
    else
        return 10 * cos(((idx-1) / 2) * 2 * π / size) + 5 - 5 * sin(((idx-1) / 2) * 2 * π / size) * im
    end
end

spec2 = zeros(ComplexF64, size)
Spectrum!(spec2, f2, size)


Am2 = spzeros(Float64, size, size)
initMat!(Am2, diag_l, diag_u, size; scale=0.1, sparsity=0.005)
@time genMat2 = nonsym(size, diag_l, diag_u, spec2, Am2, nilp)

@info "sparsity = " nnz(genMat2) / (size * size)

genspec2 = eigen!(Matrix(genMat2)).values
p2=scatter(real(spec2), imag(spec2), markercolor = :green, markersize = 7, label="given spectrum")
p2=scatter!(real(genspec2), imag(genspec2), markercolor = :red, markersize = 5, label="spectrum of generated matrix")


## Non symmetric 2
function f3(idx::Integer, size = size)
    return 10 * cos(((idx-1) / 2) * 2 * π / size) + 5
end

spec3 = zeros(ComplexF64, size)
Spectrum!(spec3, f3, size)


Am3 = spzeros(Float64, size, size)
initMat!(Am3, diag_l, diag_u, size; scale=0.1, sparsity=0.005)
@time genMat3 = nonsym(size, diag_l, diag_u, spec3, Am3, nilp)

@info "sparsity = " nnz(genMat3) / (size * size)

genspec3 = eigen!(Matrix(genMat3)).values
p3=scatter(real(spec3), imag(spec3), ylims=(-1,1), markercolor = :green, markersize = 7, label="given spectrum")
p3=scatter!(real(genspec3), imag(genspec3), ylims=(-1,1), markercolor = :red, markersize = 5, label="spectrum of generated matrix")

p4=spy(dropzeros(genMat2), legend = :none)

l = @layout [a  b;c d]
plot(p1, p2, p3, p4, layout = l)
