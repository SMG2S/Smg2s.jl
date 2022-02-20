struct Nilpotent{Ti<:Integer}
    nbOne::Ti
    size::Ti
    degree::Ti
    nilpMat::SparseMatrixCSC
end

function factorial(num::Tv, depth::Ti) where {Tv <: Real, Ti <: Integer}
    fact = num
    for i = 1:depth-1
        fact *= (num-i)
    end

    return fact
end

"""
    Nilp(nbOne::Ti, size::Ti) where {Ti<:Integer}

Create a nilpotent matrix object with parameters `nbOne` and `size`. `size` refers to the size of
nilpotent matrix to be generated, and `nbOne` refers to the number of continuous `1` on the non-zero
diagonal of generated nilptent matrix.

This is the default type of nilptent matrix used by SMG2S, in which the non-zero diagonal is selected
as the one of offset `1`. This diagonal starts with `nbOne` number of `1`, then a `0`, then `nbOne` number
of `1`, ... until to the end.

## Examples
```jldoctest; setup = :(using SMG2S)
julia> SMG2S.Nilp(2,8)
Nilpotent{Int64}(2, 8, 3,
  ⋅   1.0   ⋅    ⋅    ⋅    ⋅    ⋅    ⋅
  ⋅    ⋅   1.0   ⋅    ⋅    ⋅    ⋅    ⋅
  ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅
  ⋅    ⋅    ⋅    ⋅   1.0   ⋅    ⋅    ⋅
  ⋅    ⋅    ⋅    ⋅    ⋅   1.0   ⋅    ⋅
  ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅
  ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅   1.0
  ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅ )
```
"""
function Nilp(nbOne::Ti, size::Ti) where {Ti<:Integer}

    degree::Ti = nbOne + 1

    if nbOne > size - 1
        error("for constructing niloptent matrix, nbOne should ≦ size - 1 ")
    end

    nilpMat = spzeros(size, size)

    for i = 1:size-1
        if (i % (nbOne + 1)) != 0
            nilpMat[i, i+1] = 1.0
        end
    end

    return Nilpotent(nbOne, size, degree, nilpMat)
end

"""
    Nilp(vector::AbstractVector, size::Ti) where {Ti<:Integer}

Create a nilpotent matrix whose dimension is `size` from user-provided vector. The
non-zero diagonal of nilpotent matrix is fixed as the one with offset `1`. Therefore, the size
of given `vector` should at least be `size-1`.

## Examples
```jldoctest; setup = :(using SMG2S)
julia> vec=[1; 1; 0; 1; 1; 1; 0]
7-element Vector{Int64}:
 1
 1
 0
 1
 1
 1
 0

julia> SMG2S.Nilp(vec, 8)
Nilpotent{Int64}(3, 8, 4,
  ⋅   1.0   ⋅    ⋅    ⋅    ⋅    ⋅    ⋅
  ⋅    ⋅   1.0   ⋅    ⋅    ⋅    ⋅    ⋅
  ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅
  ⋅    ⋅    ⋅    ⋅   1.0   ⋅    ⋅    ⋅
  ⋅    ⋅    ⋅    ⋅    ⋅   1.0   ⋅    ⋅
  ⋅    ⋅    ⋅    ⋅    ⋅    ⋅   1.0   ⋅
  ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅
  ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅ )
```
"""
function Nilp(vector::AbstractVector, size::Ti) where {Ti<:Integer}

    #check
    for i = 1:size-1
        if vector[i] != 1.0 && vector[i] != 0
            vector[i] /= vector[i]
        end
    end

    nbOne::Ti = 0
    cnt::Ti = 0

    for i = 1:size-1
        if vector[i] == 1.0
            cnt += 1
        elseif vector[i] == 0.0
            if cnt > nbOne
                nbOne = cnt
            end
            cnt = 0
        end
    end

    degree::Ti = nbOne + 1

    nilpMat = spzeros(size, size)

    for i = 1:size-1
        nilpMat[i, i+1] = vector[i]
    end

    return Nilpotent(nbOne, size, degree, dropzeros(nilpMat))

end

"""
    Nilp(matrix::SparseMatrixCSC{Tv, Ti} ,size::Ti; maxdegree::Ti=80) where{Tv <: Real, Ti <: Integer}

Create a nilpotent matrix whose dimension is `size` from user-provided sparse matrix `matrix`.
This matrix should be nilpotent, and its nilpotency should not be larger than `80`. Before generation,
SMG2S.jl will check if the user-provided matrix is nilpotent matrix with nilpotency no larger than `80`.
If it doesn't satisfy any of the two, an error message will appear.

## Examples
```jldoctest; setup = :(using SMG2S; using SparseArrays)
julia> nilpMatrix = sparse([0. 1 1 0; 0 0 1 1 ; 0 0 0 1 ; 0 0 0 0])
4×4 SparseMatrixCSC{Float64, Int64} with 5 stored entries:
  ⋅   1.0  1.0   ⋅
  ⋅    ⋅   1.0  1.0
  ⋅    ⋅    ⋅   1.0
  ⋅    ⋅    ⋅    ⋅

julia> SMG2S.Nilp(nilpMatrix,4)
┌ Info: the degree of given nilpotent matrix is:
└   degree = 4
Nilpotent{Int64}(1, 4, 4,
  ⋅   1.0  1.0   ⋅
  ⋅    ⋅   1.0  1.0
  ⋅    ⋅    ⋅   1.0
  ⋅    ⋅    ⋅    ⋅ )
```

## Examples
```julia
julia> nilpMatrix = sparse([1. 1 1 0; 0 0 1 1 ; 0 0 0 1 ; 0 0 0 0])
4×4 SparseMatrixCSC{Float64, Int64} with 6 stored entries:
 1.0  1.0  1.0   ⋅
  ⋅    ⋅   1.0  1.0
  ⋅    ⋅    ⋅   1.0
  ⋅    ⋅    ⋅    ⋅

julia> SMG2S.Nilp(nilpMatrix,4)
ERROR: the given nilpotent matrix is invalid
```

"""
function Nilp(matrix::SparseMatrixCSC{Tv, Ti} ,size::Ti; maxdegree::Ti=80) where{Tv <: Real, Ti <: Integer}

    nbOne::Ti = 1
    degree::Ti = 1

    nilpMat = matrix

    while degree <= maxdegree
        if matrix == spzeros(Tv, size, size)
            break
        end
        matrix *= nilpMat
        degree += 1
    end

    if degree > maxdegree
        error("the given nilpotent matrix is invalid")
    end

    @info "the degree of given nilpotent matrix is: " degree

    return Nilpotent(nbOne, size, degree, dropzeros(nilpMat))

end

"""
    Nilp(vec::AbstractVector, diag::Ti, size::Ti) where{Ti <: Integer}

Create a nilpotent matrix whose dimension is `size` from user-provided vector. The
non-zero diagonal of nilpotent matrix is set as the one with offset `diag`. Therefore, the size
of given `vector` should at least be `size-diag`.

## Examples
```jldoctest; setup = :(using SMG2S; using SparseArrays)
julia> vec=[1; 1; 0; 1; 1; 1; 0]
7-element Vector{Int64}:
 1
 1
 0
 1
 1
 1
 0

julia> SMG2S.Nilp(vec, 3, 8)
┌ Info: the degree of given nilpotent matrix is:
└   degree = 3
Nilpotent{Int64}(1, 8, 3,
 ⋅  ⋅  ⋅  1  ⋅  ⋅  ⋅  ⋅
 ⋅  ⋅  ⋅  ⋅  1  ⋅  ⋅  ⋅
 ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅
 ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  1  ⋅
 ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  1
 ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅
 ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅
 ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅)

```
"""
function Nilp(vec::AbstractVector, diag::Ti, size::Ti) where{Ti <: Integer}
    length = size-diag

    nilpvec = vec[1:length]

    nilpMatrix = sparse(diagm(diag=>nilpvec[1:length]))

    #@info nilpvec
    return Nilp(nilpMatrix, size)
end

"""
    Nilp(nbOne::Ti, diag::Ti, size::Ti) where{Ti <: Integer}

Create a nilpotent matrix object with parameters `nbOne`, `diag` and `size`. `size` refers to the size of
nilpotent matrix to be generated, and `nbOne` refers to the maximum number of continuous `1` on the non-zero
diagonal of generated nilptent matrix. The offset of non-zero diagonal is determined by
the parameter `diag`.

The non-zero diagonal is generated as:
1. randomly generating an interger number in the interval `[1, nbOne]`

2. this number determines the number of continous `1` for this time of sampling

3. then generating another within the same interval, which determines the number of continuous `0`

4. then re-samling for a new sequence of continuous `1`

5. ...

## Examples
```julia
julia> SMG2S.Nilp(4, 3, 15)
┌ Info: the degree of given nilpotent matrix is:
└   degree = 3
Nilpotent{Int64}(1, 15, 3,
 ⋅  ⋅  ⋅  1  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅
 ⋅  ⋅  ⋅  ⋅  1  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅
 ⋅  ⋅  ⋅  ⋅  ⋅  1  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅
 ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  1  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅
 ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅
 ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅
 ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅
 ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅
 ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  1  ⋅  ⋅  ⋅
 ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  1  ⋅  ⋅
 ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  1  ⋅
 ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅
 ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅
 ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅
 ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅)
```
"""
function Nilp(nbOne::Ti, diag::Ti, size::Ti) where{Ti <: Integer}
    length = size-diag

    nilpvec = zeros(Ti, length)

    cnt = 0

    while cnt < length
        l₁ = rand(1:nbOne)
        l₀ = rand(1:nbOne)

        for i = cnt+1:l₁+cnt
            if i <= length
                nilpvec[i] = 1
            end
        end

        cnt += l₀ + l₁

    end

    return Nilp(nilpvec, diag, size)
end
