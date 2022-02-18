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

function Nilp(matrix::SparseMatrixCSC{Tv, Ti} ,size::Ti; maxdegree::Ti=80) where{Tv <: Real, Ti <: Integer}

    diagPosition::Ti = 2
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

    return Nilpotent(diagPosition, nbOne, size, degree, dropzeros(nilpMat))

end

function Nilp(vec::AbstractVector, diag::Ti, size::Ti) where{Ti <: Integer}
    length = size-diag

    nilpvec = vec[1:length]

    nilpMatrix = sparse(diagm(diag=>nilpvec[1:length]))

    #@info nilpvec
    return Nilp(nilpMatrix, size)
end

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
