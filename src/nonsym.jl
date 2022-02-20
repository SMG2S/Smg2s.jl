"""
    nonsym(nbOne::Ti, size::Ti, diag_l::Ti, diag_u::Ti, spectrum::AbstractVector{Tv}) where {Tv<:Complex, Ti<:Integer}

Generate a non symmetric matrix with default nilpotent matrix and default initialization of matrix.

The usage of `nonsym` is the same as the one of `nonherm`, please refers to that part for more detais.

Please keep in mind the constraints of given spectrum vector and initialization of matrix, click [here](getting_started.md) for more details.

"""
function nonsym(
    nbOne::Ti,
    size::Ti,
    diag_l::Ti,
    diag_u::Ti,
    spectrum::AbstractVector{Tv}
) where {Tv<:Complex, Ti<:Integer}

    Tv2 = eltype(real(spectrum[1]))
    Am = spzeros(Tv2, size, size)
    initMat!(Am, diag_l, diag_u, size)
    nilp = Nilp(nbOne, size)

    return nonsym(size, diag_l, diag_u, spectrum, Am, nilp)

end

"""
    nonsym(nbOne::Ti, size::Ti, diag_l::Ti, diag_u::Ti, spectrum::AbstractVector{Tv}, Am::SparseMatrixCSC{Tv2, Ti},) where {Tv<:Complex, Tv2<:Real, Ti<:Integer}

Generate a non symmetric matrix with default nilpotent matrix and user-provided initialization of matrix.

The usage of `nonsym` is the same as the one of `nonherm`, please refers to that part for more detais.

Please keep in mind the constraints of given spectrum vector and initialization of matrix, click [here](getting_started.md) for more details.

"""
function nonsym(
    nbOne::Ti,
    size::Ti,
    diag_l::Ti,
    diag_u::Ti,
    spectrum::AbstractVector{Tv},
    Am::SparseMatrixCSC{Tv2, Ti},
) where {Tv<:Complex, Tv2<:Real, Ti<:Integer}

    nilp = Nilp(nbOne, size)

    return nonsym(size, diag_l, diag_u, spectrum, Am, nilp)

end

"""
    nonsym(size::Ti, diag_l::Ti, diag_u::Ti, spectrum::AbstractVector{Tv}, Am::SparseMatrixCSC{Tv2, Ti}, nilp::Nilpotent{Ti}) where {Tv<:Complex, Tv2<:Real, Ti<:Integer}

Generate a non symmetric matrix with user-provided nilpotent matrix and initialization of matrix.

The usage of `nonsym` is the same as the one of `nonherm`, please refers to that part for more detais.

Please keep in mind the constraints of given spectrum vector and initialization of matrix, click [here](getting_started.md) for more details.

"""
function nonsym(
    size::Ti,
    diag_l::Ti,
    diag_u::Ti,
    spectrum::AbstractVector{Tv},
    Am::SparseMatrixCSC{Tv2, Ti},
    nilp::Nilpotent{Ti}
) where {Tv<:Complex, Tv2<:Real, Ti<:Integer}

    checkSpectrum(spectrum, size)

    for i = 1:size
        Am[i, i] = real(spectrum[i])
    end

    idx = 1
    while idx < size
        if imag(spectrum[idx]) == 0
            step = 1
        else
            step = 2
        end
        Am[idx, idx+1] = abs(imag(spectrum[idx]))
        Am[idx+1,idx] = -abs(imag(spectrum[idx]))
        idx += step
    end

    matAop = Am

    fact = 1
    fact = factorial(Tv2(2 * (nilp.degree-1)),(2 * (nilp.degree-1)))
    Am = fact * Am

    for k = 1 : 2 * (nilp.degree-1)
        matAop = matAop * nilp.nilpMat - nilp.nilpMat * matAop
        fact /= factorial(Tv2(k + 1),k+1)
        Am += fact * matAop
    end

    fact = factorial(Tv2(2 * (nilp.degree-1)), 2 * (nilp.degree-1))
    Am = (inv(fact)) * Am

    return Am
end
