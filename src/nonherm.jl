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
