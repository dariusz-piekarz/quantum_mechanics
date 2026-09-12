using QuantumMechanics

using SparseArrays


struct Grid{T<:Real}
    x::AbstractRange{T}
    y::AbstractRange{T}
    z::AbstractRange{T}

    function Grid{T}(x::AbstractRange{T}, y::AbstractRange{T}, z::AbstractRange{T}) where {T <: Real}
        if length(x) < 1 || length(y) < 1 || length(z) < 1
            throw(DimensionMismatch("Grid dimensions must be positive."))
        end

        if length(x) != length(y) || length(y) != length(z)
            throw(DimensionMismatch("Grid must have the same number of points in every dimension."))
        end

        steps = (step(x), step(y), step(z))
        if !(steps[1] == steps[2] == steps[3])
            throw(DimensionMismatch("Grid must use the same spatial step in every dimension to match the Laplacian discretization."))
        end

        return new(x, y, z)
    end
end


Grid(x::AbstractRange{T}, y::AbstractRange{T}, z::AbstractRange{T}) where {T <: Real} = Grid{T}(x, y, z)


function Grid(n::Integer, h::T) where {T <: Real}
    n < 1 && throw(DimensionMismatch("n cannot be less than 1!"))
    axis = range(zero(T), step = h, length = n)
    return Grid(axis, axis, axis)
end


function Grid(n::Integer, h::T, ::Type{T}) where {T <: Real}
    return Grid(n, h)
end


function grid_step(grid::Grid{T})::T where {T <: Real}
    return step(grid.x)
end


function grid_size(grid::Grid)
    return (length(grid.x), length(grid.y), length(grid.z))
end


function kinetic_operator(
    n::Integer,
    cut_type::Val,
    h::T,
    ::Type{T})::SparseMatrixCSC{T} where T <: Real
    return -one(T) / T(2) * laplacian(n, cut_type, h, T, 3)
end


function kinetic_operator(
    grid::Grid{T},
    cut_type::Val,
    ::Type{T})::SparseMatrixCSC{T} where T <: Real
    return -one(T) / T(2) * laplacian(length(grid.x), cut_type, grid_step(grid), T, 3)
end


function potential_operator(
    V::AbstractPotential
)(grid::Grid{T}, ::Type{T})::SparseMatrixCSC{T} where T <: Real

    diag = Vector{T}(undef, length(grid.x) * length(grid.y) * length(grid.z))
    k = 1

    for x in grid.x, y in grid.y, z in grid.z
        diag[k] = V(x, y, z)
        k += 1
    end

    return spdiagm(0 => diag)
end

