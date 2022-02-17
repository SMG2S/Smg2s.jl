function Spectrum!(spec::AbstractVector{Tv}, f::Function, size::Ti) where {Tv<:Complex, Ti<:Integer}
    for i = 1:size
        spec[i] = Tv(f(i))
    end
end

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
