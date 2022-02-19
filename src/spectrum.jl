"""
    Spectrum!(spec::AbstractVector{Tv}, f::Function, size::Ti) where {Tv<:Complex, Ti<:Integer}

Set the spectrum from user provided function.

## Examples
```jldoctest; setup = :(using SMG2S; using LinearAlgebra)
julia> function f(idx::Integer, size = 10)
           return 10 * cos(((idx-1) / 2) * 2 * π / size) + sin(((idx-1) / 2) * 2 * π / size) * im
       end
f (generic function with 2 methods)

julia> vec=zeros(ComplexF64, 10);

julia> Spectrum!(vec, f, 10);

julia> vec
10-element Vector{ComplexF64}:
                  10.0 + 0.0im
     9.510565162951535 + 0.3090169943749474im
     8.090169943749475 + 0.5877852522924731im
     5.877852522924732 + 0.8090169943749475im
    3.0901699437494745 + 0.9510565162951535im
 6.123233995736766e-16 + 1.0im
   -3.0901699437494736 + 0.9510565162951536im
     -5.87785252292473 + 0.8090169943749475im
    -8.090169943749473 + 0.5877852522924732im
    -9.510565162951535 + 0.3090169943749475im

```

"""
function Spectrum!(spec::AbstractVector{Tv}, f::Function, size::Ti) where {Tv<:Complex, Ti<:Integer}
    for i = 1:size
        spec[i] = Tv(f(i))
    end
end

"""
    Spectrum!(spec::AbstractVector{Tv}, vec::AbstractVector{Tv}, size::Ti) where {Tv<:Complex, Ti<:Integer}

Set the spectrum from user provided vector.

"""
function Spectrum!(spec::AbstractVector{Tv}, vec::AbstractVector{Tv}, size::Ti) where {Tv<:Complex, Ti<:Integer}
    for i = 1:size
        spec[i] = Tv(vec[i])
    end
end

##only for non-symmetric cases
function checkSpectrum(spectrum::AbstractVector{Tv}, size::Ti) where {Tv <: Complex, Ti <: Integer}
    idx = 1

    while idx < size
        if imag(spectrum[idx]) == 0
            step = 1
        else
            if (imag(spectrum[idx]) != -imag(spectrum[idx+1])) && (real(spectrum[idx]) != real(spectrum[idx+1]))
                @info spectrum[idx], spectrum[idx+1],conj(spectrum[idx+1])
                error(
                    "for initialisationg of non-symmetric matrix, the given spectrum is invalid, please follow the instruction",
                    )
            end
            step = 2
        end

        idx += step
    end
end
