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
