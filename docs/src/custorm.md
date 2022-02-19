# Customization

## Nilpotent Matrix

### Default Nilpotent Matrix

```@docs
Nilp(nbOne::Ti, size::Ti) where {Ti<:Integer}
```

### Construction from User-provided vector

```@docs
Nilp(vector::AbstractVector, size::Ti) where {Ti<:Integer}
```

### Construction from User-provided Sparse Matrix

```@docs
Nilp(matrix::SparseMatrixCSC{Tv, Ti} ,size::Ti; maxdegree::Ti=80) where{Tv <: Real, Ti <: Integer}
```

### Construction from User-provided vector with specific diagonal offset

```@docs
Nilp(vec::AbstractVector, diag::Ti, size::Ti) where{Ti <: Integer}
```

### Construction from specific diagonal offset

```@docs
Nilp(nbOne::Ti, diag::Ti, size::Ti) where{Ti <: Integer}
```

## Initialization of Matrix

### Non-Hermtian Case

```@docs
initMat!(matrix::SparseMatrixCSC{Tv, Ti}, diag_l::Ti, diag_u::Ti, size::Ti; scale::Real = 1.0, shift::Real = 0.0, sparsity::Real = 0.9) where {Tv<:Complex, Ti<:Integer}
```

### Non-Symmetric Case

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

## Assembling the Customizations

### Non-Hermitian Case

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
