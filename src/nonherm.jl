"""
    nonherm(nbOne::Ti, size::Ti, diag_l::Ti, diag_u::Ti, spectrum::AbstractVector{Tv}) where {Tv<:Complex, Ti<:Integer}

Generate a non hermitian matrix with default nilpotent matrix and default initialization of matrix.

## Examples
```julia
julia> spec=[1+im, 2, 4, 6, 7-2im, 4+3im, 7, 8, 9, 10.2-im]
10-element Vector{ComplexF64}:
  1.0 + 1.0im
  2.0 + 0.0im
  4.0 + 0.0im
  6.0 + 0.0im
  7.0 - 2.0im
  4.0 + 3.0im
  7.0 + 0.0im
  8.0 + 0.0im
  9.0 + 0.0im
 10.2 - 1.0im

julia> nonherm(3, 10, -8, -2, spec)
10×10 SparseMatrixCSC{ComplexF64, Int64} with 66 stored entries:
   1.04657+1.0im    -0.4959+0.5im  0.0833701+0.0833333im  0.00347223+0.00347222im            ⋅                ⋅               ⋅                     ⋅                ⋅           ⋅
 -0.231523+0.0im     1.9567+0.0im   -1.00012+0.0im        3.02285e-5+0.0im                   ⋅                ⋅               ⋅                     ⋅                ⋅           ⋅
  0.207446+0.0im  -0.115496+0.0im    3.94688+0.0im          -1.00394+0.0im                   ⋅                ⋅               ⋅                     ⋅                ⋅           ⋅
  0.766797+0.0im   0.949636+0.0im   0.347018+0.0im           6.04985+0.0im                   ⋅                ⋅               ⋅                     ⋅                ⋅           ⋅
  0.413475+0.0im   0.847118+0.0im    1.11204+0.0im          0.240384+0.0im           7.05096-2.0im    1.50334-2.5im  0.500002-0.666667im  0.0277777-0.0381944im      ⋅           ⋅
 -0.114481+0.0im   0.314934+0.0im   0.463301+0.0im          0.257056+0.0im         -0.253119+0.0im     3.9802+3.0im  -1.49669+1.5im       -0.166576+0.25im           ⋅           ⋅
  0.592639+0.0im   0.167683+0.0im  -0.318934+0.0im          0.847507+0.0im          0.225334+0.0im  -0.292356+0.0im    6.8867+0.0im       -0.507753+0.0im            ⋅           ⋅
  0.112486+0.0im   0.332264+0.0im   0.799049+0.0im          0.349224+0.0im          0.842714+0.0im    1.37186+0.0im  0.545475+0.0im         8.08213+0.0im            ⋅           ⋅
   0.96467+0.0im    1.36421+0.0im    1.13479+0.0im          0.769023+0.0im         -0.219469+0.0im  -0.196588+0.0im   0.26921+0.0im        0.154841+0.0im        9.0+0.0im  -0.6+0.5im
           ⋅        0.12781+0.0im   0.438092+0.0im           1.12886+0.0im          0.438939+0.0im   0.466332+0.0im  0.160009+0.0im        0.022096+0.0im            ⋅      10.2-1.0im
```

"""
function nonherm(
    nbOne::Ti,
    size::Ti,
    diag_l::Ti,
    diag_u::Ti,
    spectrum::AbstractVector{Tv}
) where {Tv<:Complex, Ti<:Integer}

    Am = spzeros(Tv, size, size)
    initMat!(Am, diag_l, diag_u, size)

    nilp = Nilp(nbOne, size)
    return nonherm(size, diag_l, diag_u, spectrum, Am, nilp)

end

"""
    nonherm(nbOne::Ti, size::Ti, diag_l::Ti, diag_u::Ti, spectrum::AbstractVector{Tv}, Am::SparseMatrixCSC{Tv, Ti}) where {Tv<:Complex, Ti<:Integer}

Generate a non hermitian matrix with default nilpotent matrix and user-provided initialization of matrix.

## Examples
```julia
julia> spec=[1+im, 2, 4, 6, 7-2im, 4+3im, 7, 8, 9, 10.2-im]
10-element Vector{ComplexF64}:
  1.0 + 1.0im
  2.0 + 0.0im
  4.0 + 0.0im
  6.0 + 0.0im
  7.0 - 2.0im
  4.0 + 3.0im
  7.0 + 0.0im
  8.0 + 0.0im
  9.0 + 0.0im
 10.2 - 1.0im

julia> Am = spzeros(ComplexF64, 10, 10);
julia> initMat!(Am, -8, -2, 10; scale=0.1, sparsity=0.005)
julia> nonherm(3, 10, -8, -2, spec, Am)
10×10 SparseMatrixCSC{ComplexF64, Int64} with 22 stored entries:
 1.0+1.0im  -0.5+0.5im  0.0833333+0.0833333im  0.00347222+0.00347222im      ⋅          ⋅           ⋅                     ⋅                ⋅           ⋅
     ⋅       2.0+0.0im       -1.0+0.0im                   ⋅                 ⋅          ⋅           ⋅                     ⋅                ⋅           ⋅
     ⋅           ⋅            4.0+0.0im              -1.0+0.0im             ⋅          ⋅           ⋅                     ⋅                ⋅           ⋅
     ⋅           ⋅                ⋅                   6.0+0.0im             ⋅          ⋅           ⋅                     ⋅                ⋅           ⋅
     ⋅           ⋅                ⋅                       ⋅             7.0-2.0im  1.5-2.5im   0.5-0.666667im  0.0277778-0.0381944im      ⋅           ⋅
     ⋅           ⋅                ⋅                       ⋅                 ⋅      4.0+3.0im  -1.5+1.5im       -0.166667+0.25im           ⋅           ⋅
     ⋅           ⋅                ⋅                       ⋅                 ⋅          ⋅       7.0+0.0im            -0.5+0.0im            ⋅           ⋅
     ⋅           ⋅                ⋅                       ⋅                 ⋅          ⋅           ⋅                 8.0+0.0im            ⋅           ⋅
     ⋅           ⋅                ⋅                       ⋅                 ⋅          ⋅           ⋅                     ⋅            9.0+0.0im  -0.6+0.5im
     ⋅           ⋅                ⋅                       ⋅                 ⋅          ⋅           ⋅                     ⋅                ⋅      10.2-1.0im
```

"""
function nonherm(
    nbOne::Ti,
    size::Ti,
    diag_l::Ti,
    diag_u::Ti,
    spectrum::AbstractVector{Tv},
    Am::SparseMatrixCSC{Tv, Ti}
) where {Tv<:Complex, Ti<:Integer}

    Tv2 = eltype(real(spectrum[1]))

    nilp = Nilp(nbOne, size)

    return nonherm(size, diag_l, diag_u, spectrum, Am, nilp)

end

"""
    nonherm(size::Ti, diag_l::Ti, diag_u::Ti, spectrum::AbstractVector{Tv}, Am::SparseMatrixCSC{Tv, Ti}, nilp::Nilpotent{Ti}) where {Tv<:Complex, Ti<:Integer}

Generate a non hermitian matrix with user-provided nilpotent matrix and initialization of matrix.

## Examples
```julia
julia> spec=[1+im, 2, 4, 6, 7-2im, 4+3im, 7, 8, 9, 10.2-im]
10-element Vector{ComplexF64}:
  1.0 + 1.0im
  2.0 + 0.0im
  4.0 + 0.0im
  6.0 + 0.0im
  7.0 - 2.0im
  4.0 + 3.0im
  7.0 + 0.0im
  8.0 + 0.0im
  9.0 + 0.0im
 10.2 - 1.0im

julia> Am = spzeros(ComplexF64, 10, 10);
julia> initMat!(Am, -8, -2, 10; scale=0.1, sparsity=0.005)
julia> nilp=SMG2S.Nilp(4, 3, 10)
┌ Info: the degree of given nilpotent matrix is:
└   degree = 2
Nilpotent{Int64}(1, 10, 2,
 ⋅  ⋅  ⋅  1  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅
 ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅
 ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅
 ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅
 ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  1  ⋅  ⋅
 ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅
 ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  1
 ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅
 ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅
 ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅)

julia> nonherm(10, -8, -2, spec, Am, nilp)
10×10 SparseMatrixCSC{ComplexF64, Int64} with 13 stored entries:
 1.0+1.0im      ⋅          ⋅      -2.5+0.5im      ⋅          ⋅          ⋅           ⋅          ⋅           ⋅
     ⋅      2.0+0.0im      ⋅           ⋅          ⋅          ⋅          ⋅           ⋅          ⋅           ⋅
     ⋅          ⋅      4.0+0.0im       ⋅          ⋅          ⋅          ⋅           ⋅          ⋅           ⋅
     ⋅          ⋅          ⋅       6.0+0.0im      ⋅          ⋅          ⋅           ⋅          ⋅           ⋅
     ⋅          ⋅          ⋅           ⋅      7.0-2.0im      ⋅          ⋅      -0.5-1.0im      ⋅           ⋅
     ⋅          ⋅          ⋅           ⋅          ⋅      4.0+3.0im      ⋅           ⋅          ⋅           ⋅
     ⋅          ⋅          ⋅           ⋅          ⋅          ⋅      7.0+0.0im       ⋅          ⋅      -1.6+0.5im
     ⋅          ⋅          ⋅           ⋅          ⋅          ⋅          ⋅       8.0+0.0im      ⋅           ⋅
     ⋅          ⋅          ⋅           ⋅          ⋅          ⋅          ⋅           ⋅      9.0+0.0im       ⋅
     ⋅          ⋅          ⋅           ⋅          ⋅          ⋅          ⋅           ⋅          ⋅      10.2-1.0im
```
"""
function nonherm(
    size::Ti,
    diag_l::Ti,
    diag_u::Ti,
    spectrum::AbstractVector{Tv},
    Am::SparseMatrixCSC{Tv, Ti},
    nilp::Nilpotent{Ti}
) where {Tv<:Complex, Ti<:Integer}

    Tv2 = eltype(real(spectrum[1]))
    nbOne = nilp.nbOne
    for i = 1:size
        Am[i, i] = spectrum[i]
    end

    matAop = Am

    fact = 1
    fact = factorial(Tv2(2 * (nilp.degree-1)), 2 * (nilp.degree-1))
    Am = fact * Am

    for k = 1 : 2 * (nilp.degree-1)
        matAop = matAop * nilp.nilpMat - nilp.nilpMat * matAop
        fact /= factorial(Tv2(k + 1),k+1)
        Am += fact * matAop
    end

    fact = factorial(Tv2(2 * (nilp.degree-1)),2 * (nilp.degree-1))
    Am = (inv(fact)) * Am

    return Am
end
