
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
