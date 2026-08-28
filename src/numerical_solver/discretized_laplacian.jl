using SparseArrays


function raw_laplacian(n::Integer, h::T, ::Type{T})::SparseMatrixCSC{T} where T <: Real
    n <= 0 && throw(DimensionMismatch("n cannot be less than 1!"))
    
    rows = repeat(1:n, inner = [3])
    cols = rows + repeat(0:2, outer=[n])
    values = repeat(T[1, -2, 1], outer=[n])
    return T(1 / h ^ 2) * sparse(rows, cols, values)
end


function laplacian(
    n::Integer,
    ::Val{:start},
    h::T,
    ::Type{T}
) where T <: Real
    raw = raw_laplacian(n, h, T)
    return raw[:, 3:(n + 2)]
end


function laplacian(
    n::Integer,
    ::Val{:end},
    h::T,
    ::Type{T}
) where T <: Real

    raw = raw_laplacian(n, h, T)
    return raw[:, 1:n]
end


function laplacian(
    n::Integer,
    ::Val{:twotail},
    h::T,
    ::Type{T}
) where T <: Real

    raw = raw_laplacian(n, h, T)
    return raw[:, 2:(n + 1)]
end