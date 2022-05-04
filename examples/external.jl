using Smg2s
using Test
using SparseArrays
using Plots

using LinearAlgebra
using Random
using DelimitedFiles

V1 = readdlm("./examples/data/v1.txt", skipstart=2)[:,2]
V2 = readdlm("./examples/data/v2.txt", skipstart=2)[:,2:3]

@test size(V1,1) == size(V2,1)

n = size(V1,1)

function f1(idx::Integer, n = n)
    return V1[idx] + 0.0 * im
end

function f2(idx::Integer, n = n)
    return V2[idx][1] + V2[idx][1] * im
end

##initalization of SMG2S, which effects only the sparsity of generated matrices
diag_l = -50
diag_u = -10

nbOne = 6

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

#nilpMatrix = sparse(diagm(10=>nilpvec[1:n-10]))

nilp = Nilp(nilpvec, 10, n)

##v1.txt case: non-symmetric matrix
spec1 = zeros(ComplexF64, n)
Spectrum!(spec1, f1, n)

Am = spzeros(Float64, n, n)
initMat!(Am, diag_l, diag_u, n; scale=0.1, sparsity=.5)
@time genMat1 = nonsym(n, diag_l, diag_u, spec1, Am, nilp)

@info "sparsity = " nnz(genMat1) / (n * n)

genspec1 = eigen!(Matrix(genMat1)).values
p1=scatter(real(spec1), imag(spec1), ylims=(-1,1), markercolor = :green, markern = 7, markersize = 4, label="given spectrum")
p1=scatter!(real(genspec1), imag(genspec1), ylims=(-1,1), markercolor = :red, markern = 5, markersize = 3, label="spectrum of generated matrix",legendfontsize=6,legend=:topleft)

##v2.txt case: non-Hermitian matrix
spec2 = zeros(ComplexF64, n)
Spectrum!(spec2, f2, n)

Am = spzeros(ComplexF64, n, n)
initMat!(Am, diag_l, diag_u, n; scale=0.1, sparsity=.005)
@time genMat2 = nonherm(n, diag_l, diag_u, spec2, Am, nilp)

@info "sparsity = " nnz(genMat2) / (n * n)

genspec2 = eigen!(Matrix(genMat2)).values
p2=scatter(real(spec2), imag(spec2), markercolor = :green, markern = 7, markersize = 4, label="given spectrum")
p2=scatter!(real(genspec2), imag(genspec2), markercolor = :red, markern = 5, markersize = 3, label="spectrum of generated matrix", legendfontsize=5, legend=:topleft)

p3=spy(dropzeros(genMat1), legend = :none)

l = @layout [a  b;c]
plot(p1, p2, p3, layout = l)
