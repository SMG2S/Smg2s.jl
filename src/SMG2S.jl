module SMG2S

using LinearAlgebra
using SparseArrays
using Random

export Nilpotent, Nilp
export Spectrum!, checkSpectrum
export initMat!
export nonherm, nonsym

include("nilpotent.jl")
include("spectrum.jl")
include("init.jl")
include("nonsym.jl")
include("nonherm.jl")

#=
struct Nilpotent{Ti<:Integer}
    diagPosition::Ti
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

function Nilp(nbOne::Ti, size::Ti) where {Ti<:Integer}

    diagPosition::Ti = 2
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

    return Nilpotent(diagPosition, nbOne, size, degree, nilpMat)
end

function Nilp(vector::AbstractVector, size::Ti) where {Ti<:Integer}
    diagPosition::Ti = 2

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

    return Nilpotent(diagPosition, nbOne, size, degree, dropzeros(nilpMat))

end

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

function checkSpectrum(spectrum::AbstractVector{Tv}, size::Ti) where {Tv <: Complex, Ti <: Integer}
    idx = 1

    while idx < size
        if imag(spectrum[idx]) == 0
            step = 1
        else
            if (imag(spectrum[idx]) != imag(spectrum[idx+1])) && (real(spectrum[idx]) != real(spectrum[idx+1]))
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

function nonherm(
    nbOne::Ti,
    size::Ti,
    diag_l::Ti,
    diag_u::Ti,
    spectrum::AbstractVector{Tv}
) where {Tv<:Complex, Ti<:Integer}

    Tv2 = eltype(real(spectrum[1]))
    Am = spzeros(Tv, size, size)
    initMat!(Am, diag_l, diag_u, size)

    nilp = Nilp(nbOne, size)
    for i = 1:size
        Am[i, i] = spectrum[i]
    end

    matAop = Am

    fact = 1
    fact = factorial(Tv2(2 * nilp.nbOne),2 * nilp.nbOne)
    Am = fact * Am

    for k = 1 : 2 * nilp.nbOne
        matAop = matAop * nilp.nilpMat - nilp.nilpMat * matAop
        fact /= factorial(Tv2(k + 1), k+1)
        Am += fact * matAop
    end

    fact = factorial(Tv2(2 * nilp.nbOne), 2 * nilp.nbOne)
    Am = (inv(fact)) * Am

    return Am
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
    for i = 1:size
        Am[i, i] = spectrum[i]
    end

    matAop = Am

    fact = 1
    fact = factorial(Tv2(2 * nilp.nbOne),2 * nilp.nbOne)
    Am = fact * Am

    for k = 1 : 2 * nilp.nbOne
        matAop = matAop * nilp.nilpMat - nilp.nilpMat * matAop
        fact /= factorial(Tv2(k + 1),k+1)
        Am += fact * matAop
    end

    fact = factorial(Tv2(2 * nilp.nbOne),2 * nilp.nbOne)
    Am = (inv(fact)) * Am

    return Am
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
    fact = factorial(Tv2(2 * nilp.nbOne), 2 * nilp.nbOne)
    Am = fact * Am

    for k = 1 : 2 * nilp.nbOne
        matAop = matAop * nilp.nilpMat - nilp.nilpMat * matAop
        fact /= factorial(Tv2(k + 1),k+1)
        Am += fact * matAop
    end

    fact = factorial(Tv2(2 * nilp.nbOne),2 * nilp.nbOne)
    Am = (inv(fact)) * Am

    return Am
end


function nonsym(
    nbOne::Ti,
    size::Ti,
    diag_l::Ti,
    diag_u::Ti,
    spectrum::AbstractVector{Tv}
) where {Tv<:Complex, Ti<:Integer}

    checkSpectrum(spectrum, size)

    Tv2 = eltype(real(spectrum[1]))

    Am = spzeros(Tv2, size, size)
    initMat!(Am, diag_l, diag_u, size)

    nilp = Nilp(nbOne, size)

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
    fact = factorial(Tv2(2 * nilp.nbOne),2 * nilp.nbOne)
    Am = fact * Am

    for k = 1 : 2 * nilp.nbOne
        matAop = matAop * nilp.nilpMat - nilp.nilpMat * matAop
        fact /= factorial(Tv2(k + 1),k+1)
        Am += fact * matAop
    end

    fact = factorial(Tv2(2 * nilp.nbOne),2 * nilp.nbOne)
    Am = (inv(fact)) * Am

    return Am
end


function nonsym(
    nbOne::Ti,
    size::Ti,
    diag_l::Ti,
    diag_u::Ti,
    spectrum::AbstractVector{Tv},
    Am::SparseMatrixCSC{Tv2, Ti},
) where {Tv<:Complex, Tv2<:Real, Ti<:Integer}

    checkSpectrum(spectrum, size)

    nilp = Nilp(nbOne, size)

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
    fact = factorial(Tv2(2 * nilp.nbOne),2 * nilp.nbOne)
    Am = fact * Am

    for k = 1 : 2 * nilp.nbOne
        matAop = matAop * nilp.nilpMat - nilp.nilpMat * matAop
        fact /= factorial(Tv2(k + 1),k+1)
        Am += fact * matAop
    end

    fact = factorial(Tv2(2 * nilp.nbOne),2 * nilp.nbOne)
    Am = (inv(fact)) * Am

    return Am
end


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
    fact = factorial(Tv2(2 * nilp.nbOne),(2 * nilp.nbOne))
    Am = fact * Am

    for k = 1 : 2 * nilp.nbOne
        matAop = matAop * nilp.nilpMat - nilp.nilpMat * matAop
        fact /= factorial(Tv2(k + 1),k+1)
        Am += fact * matAop
    end

    fact = factorial(Tv2(2 * nilp.nbOne), 2 * nilp.nbOne)
    Am = (inv(fact)) * Am

    return Am
end

=#
end #endmodule
