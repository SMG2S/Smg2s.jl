"""
    initMat!(matrix::SparseMatrixCSC{Tv, Ti}, diag_l::Ti, diag_u::Ti, size::Ti; scale::Real = 1.0, shift::Real = 0.0, sparsity::Real = 0.9) where {Tv<:Complex, Ti<:Integer}

Initialization of matrix to be generated. This function is for non-Hermitian case, in which the entries of matrix are complex scalars.

In this function, the parameters `diag_l` and `diag_u` determine the range (`[diag_l, diag_u]`) of lower triangular part to be filled. The two parameters should also be negatif, which
refer to the offfsets of diagonals in the lower triangular part. The parameter `size` defines the size of matrix to be generated.

The non-zero elements of the initialized matrices are randomly generated in the interval (0,1). Two optionals parameters `shift` and `scale` allow shifting and scaling this interval.

Another optional parameter `sparsity` is more important. It determines the possility of the elements within the band of diagonal determined by `diag_l` and `diag_u` to be non-zeros.
In other words, the initalized matrix would be more sparse with a smaller number of `sparisty` parameter.

## Examples
```julia
julia> Am=spzeros(ComplexF64, 50, 50);
julia> initMat!(Am, -20, -10, 50)
julia> Am
50×50 SparseMatrixCSC{ComplexF64, Int64} with 335 stored entries:
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⣄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⣟⡳⣄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⣿⣿⣢⣵⣄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠙⢟⣾⣿⣽⢲⣄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠙⢽⣿⣛⣿⣵⣄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠙⠽⣻⡾⣿⣷⣄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠙⢟⣻⢿⢻⣗⣄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠙⢿⣽⣯⣻⣳⣄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢿⣛⢿⡾⣷⣄⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠑⢯⡽⣯⣜⣷⢄⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⠚⠋⠛⠚⠓⠀⠀⠀⠀⠀

```

## Examples
```julia
julia> Am=spzeros(ComplexF64, 50, 50);
julia> initMat!(Am, -20, -10, 50, sparsity=0.1)
julia> Am
50×50 SparseMatrixCSC{ComplexF64, Int64} with 36 stored entries:
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠠⢀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠉⠂⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠂⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⢂⠀⠈⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠂⠈⠀⠤⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⢀⢀⠅⠂⠄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠂⠀⠠⠈⠰⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠂⠈⠀⠡⢀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠄⢄⡆⢀⠀⠀⠀⠀⠀⠀

```⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀

"""
function initMat!(
    matrix::SparseMatrixCSC{Tv, Ti},
    diag_l::Ti,
    diag_u::Ti,
    size::Ti;
    scale::Real = 1.0,
    shift::Real = 0.0,
    sparsity::Real = 0.9
) where {Tv<:Complex, Ti<:Integer}

    diag_l = abs(diag_l)
    diag_u = abs(diag_u)

    if diag_l < diag_u
        error(
            "for initialisationg of matrix, please ensure abs(diag_l) < abas(diag_u) ",
            )
        end

    if diag_u >= size
        error(
            "for initialisationg of matrix, please ensure abs(diag_u) < size ",
            )
    end

    cnt::Ti = (diag_l - diag_u + 1) * (2 * size - diag_l - diag_u) / 2

    rnd = rand(MersenneTwister(1234), Float64, cnt)

    for i = 1:cnt
        if real(rnd[i]) < 1.0 - sparsity
            rnd[i] = 0.0
        end
    end

    idx = 0
    for i = diag_u + 1:size
        for j = max(1, i - diag_l) : i - diag_u
            idx += 1
            matrix[i, j] = Tv(scale * (rnd[idx]) + shift)
        end
    end
end

"""
    initMat!(matrix::SparseMatrixCSC{Tv, Ti}, diag_l::Ti, diag_u::Ti, size::Ti; scale::Real = 1.0, shift::Real = 0.0, sparsity::Real = 0.9) where {Tv<:Real, Ti<:Integer}

Initialization of matrix to be generated. This function is for non-Symmetric case, in which the entries of matrix are real scalars.

The usage of this function is quite similar as the one for non-Hermitian case, the only difference is that the entries of matrices are real scalars.

## Examples
```julia
julia> Am=spzeros(Float64, 50, 50);
julia> initMat!(Am, -20, -10, 50, sparsity=0.5)
julia> Am
50×50 SparseMatrixCSC{Float64, Int64} with 178 stored entries:
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⣄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠌⡱⣄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠍⠫⢠⡀⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠈⠏⣌⠺⣑⠰⢀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠙⣣⢑⡿⢡⢀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠅⣊⡼⠖⠷⡄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠉⢀⣫⠇⢛⡆⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠑⠟⡕⣭⠉⣱⡄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠉⢐⡚⢘⠔⣡⢄⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠪⡌⡆⣜⣆⢀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠒⠉⠂⠘⠁⠀⠀⠀⠀⠀

```

"""
function initMat!(
    matrix::SparseMatrixCSC{Tv, Ti},
    diag_l::Ti,
    diag_u::Ti,
    size::Ti;
    scale::Real = 1.0,
    shift::Real = 0.0,
    sparsity::Real = 0.9
) where {Tv<:Real, Ti<:Integer}

    diag_l = abs(diag_l)
    diag_u = abs(diag_u)

    if diag_u < 2
        error(
            "for initialisationg of non-symmetric matrix, please ensure abs(diag_u) >= 2 ",
            )
    end


    if diag_l < diag_u
        error(
            "for initialisationg of matrix, please ensure abs(diag_l) < abs(diag_u) ",
            )
    end

    if diag_u >= size
        error(
            "for initialisationg of matrix, please ensure abs(diag_u) < size ",
            )
    end

    cnt::Ti = (diag_l - diag_u + 1) * (2 * size - diag_l - diag_u) / 2

    rnd = rand(MersenneTwister(1234), Float64, cnt)

    for i = 1:cnt
        if real(rnd[i]) < 1.0 - sparsity
            rnd[i] = 0.0
        end
    end

    idx = 0
    for i = diag_u + 1:size
        for j = max(1, i - diag_l) : i - diag_u
            idx += 1
            matrix[i, j] = Tv(scale * (rnd[idx]) + shift)
        end
    end

end
