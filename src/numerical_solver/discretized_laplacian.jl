using SparseArrays


function raw_laplacian_1d(
    n::Integer,
    h::T,
    ::Type{T},
    constant::T = one(T)
)::SparseMatrixCSC{T} where {T <: Real}
    n <= 0 && throw(DimensionMismatch("n cannot be less than 1!"))
    (iszero(h) || !isfinite(h)) && throw(ArgumentError("Grid spacing h must be finite and nonzero."))

    rows = repeat(1:n, inner = [3])
    cols = rows + repeat(0:2, outer = [n])
    values = repeat(T[1, -2, 1], outer = [n])
    return (constant / h^2) * sparse(rows, cols, values)
end


function laplacian_1d(
    n::Integer,
    h::T,
    ::Type{T},
    constant::T = one(T)
)::SparseMatrixCSC{T} where T <: Real

    raw = raw_laplacian_1d(n, h, T, constant)
    return raw[:, 2:(n + 1)]
end


function laplacian(
    n::Integer,
    h::T,
    ::Type{T},
    dim::Integer = 3,
    constant::T = one(T)
)::SparseMatrixCSC{T} where {T <: Real}
    dim < 1 && throw(DimensionMismatch("Cannot create the discretized Laplacian for a dimension < 1!"))

    lp = spzeros(T, n^dim, n^dim)
    lp1 = laplacian_1d(n, h, T, constant)

    for i in 1:dim
        left_id = spdiagm(0 => ones(T, n^(i - 1)))
        right_id = spdiagm(0 => ones(T, n^(dim - i)))
        lp += kron(left_id, lp1, right_id)
    end

    return lp
end
