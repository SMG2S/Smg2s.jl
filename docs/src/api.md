# API Reference

## Nilpotent Matrix

```@docs
Nilp(nbOne::Ti, size::Ti) where {Ti<:Integer}
```

```@docs
Nilp(vector::AbstractVector, size::Ti) where {Ti<:Integer}
```

```@docs
Nilp(matrix::SparseMatrixCSC{Tv, Ti} ,size::Ti; maxdegree::Ti=80) where{Tv <: Real, Ti <: Integer}
```

```@docs
Nilp(vec::AbstractVector, diag::Ti, size::Ti) where{Ti <: Integer}
```

```@docs
Nilp(nbOne::Ti, diag::Ti, size::Ti) where{Ti <: Integer}
```
## Matrix Initialization

```@docs
initMat!(matrix::SparseMatrixCSC{Tv, Ti}, diag_l::Ti, diag_u::Ti, size::Ti; scale::Real = 1.0, shift::Real = 0.0, sparsity::Real = 0.9) where {Tv<:Complex, Ti<:Integer}
```

```@docs
initMat!(matrix::SparseMatrixCSC{Tv, Ti}, diag_l::Ti, diag_u::Ti, size::Ti; scale::Real = 1.0, shift::Real = 0.0, sparsity::Real = 0.9) where {Tv<:Real, Ti<:Integer}
```


## Set Spectrum

```@docs
Spectrum!(spec::AbstractVector{Tv}, f::Function, size::Ti) where {Tv<:Complex, Ti<:Integer}
```

```@docs
Spectrum!(spec::AbstractVector{Tv}, vec::AbstractVector{Tv}, size::Ti) where {Tv<:Complex, Ti<:Integer}
```

## Generation

### Non-Hermitian Matrix

```@docs
nonherm(nbOne::Ti, size::Ti, diag_l::Ti, diag_u::Ti, spectrum::AbstractVector{Tv}) where {Tv<:Complex, Ti<:Integer}
```

```@docs
nonherm(nbOne::Ti, size::Ti, diag_l::Ti, diag_u::Ti, spectrum::AbstractVector{Tv}, Am::SparseMatrixCSC{Tv, Ti}) where {Tv<:Complex, Ti<:Integer}
```

```@docs
nonherm(size::Ti, diag_l::Ti, diag_u::Ti, spectrum::AbstractVector{Tv}, Am::SparseMatrixCSC{Tv, Ti}, nilp::Nilpotent{Ti}) where {Tv<:Complex, Ti<:Integer}
```

### Non-Symmetric Matrix

```@docs
nonsym(nbOne::Ti, size::Ti, diag_l::Ti, diag_u::Ti, spectrum::AbstractVector{Tv}) where {Tv<:Complex, Ti<:Integer}
```

```@docs
nonsym(nbOne::Ti, size::Ti, diag_l::Ti, diag_u::Ti, spectrum::AbstractVector{Tv}, Am::SparseMatrixCSC{Tv2, Ti},) where {Tv<:Complex, Tv2<:Real, Ti<:Integer}
```

```@docs
nonsym(size::Ti, diag_l::Ti, diag_u::Ti, spectrum::AbstractVector{Tv}, Am::SparseMatrixCSC{Tv2, Ti}, nilp::Nilpotent{Ti}) where {Tv<:Complex, Tv2<:Real, Ti<:Integer}
```
