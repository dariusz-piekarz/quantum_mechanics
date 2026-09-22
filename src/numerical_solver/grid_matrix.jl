using QuantumMechanics

using SparseArrays


struct GridMatrix{T<:Real}
    grid::Matrix{T}

    function GridMatrix{T}(grid::AbstractMatrix{<:Real}) where {T <: Real}
        return new{T}(Matrix{T}(grid))
    end
end


function GridMatrix(x::AbstractRange{T}, y::AbstractRange{T}, z::AbstractRange{T}) where {T <: Real}
    if length(x) < 2 || length(y) < 2 || length(z) < 2
        throw(DimensionMismatch("Grid dimensions must contain at least two points."))
    end

    if x[2] - x[1] != y[2] - y[1] || y[2] - y[1] != z[2] - z[1]
        throw(DimensionMismatch("Grid must use the same spatial step in every dimension to match the Laplacian discretization."))
    end

    grid = transpose([
                vec([i for i in x, j in y, k in z])';
                vec([j for i in x, j in y, k in z])';
                vec([k for i in x, j in y, k in z])'
                    ])
    return GridMatrix{T}(grid)
end


function GridMatrix(base_range::AbstractRange{T}) where {T <: Real}
    if length(base_range) < 2
        throw(DimensionMismatch("Grid dimensions must contain at least two points."))
    end

    grid = transpose([
                vec([i for i in base_range, j in base_range, k in base_range])';
                vec([j for i in base_range, j in base_range, k in base_range])';
                vec([k for i in base_range, j in base_range, k in base_range])'
                    ])
    return GridMatrix{T}(grid)
end


function GridMatrix(
    x_min::T,
    x_max::T,
    y_min::T,
    y_max::T,
    z_min::T,
    z_max::T,
    h::T) where {T <: Real}

    if x_max <= x_min || y_max <= y_min || z_max <= z_min
        throw(DimensionMismatch("Grid max values must be greater than min values."))
    end

    if h <= zero(T)
        throw(DimensionMismatch("Grid step size must be positive."))
    end

    x = range(x_min, step = h, length = round(Int, (x_max - x_min) / h) + 1)
    y = range(y_min, step = h, length = round(Int, (y_max - y_min) / h) + 1)
    z = range(z_min, step = h, length = round(Int, (z_max - z_min) / h) + 1)

    return GridMatrix(x, y, z)
end


function GridMatrix(
    min_v::T,
    max_v::T,
    h::T) where {T <: Real}

    if max_v <= min_v
        throw(DimensionMismatch("Grid max value must be greater than min value."))
    end

    if h <= zero(T)
        throw(DimensionMismatch("Grid step size must be positive."))
    end

    x = range(min_v, step = h, length = round(Int, (max_v - min_v) / h) + 1)
    y = range(min_v, step = h, length = round(Int, (max_v - min_v) / h) + 1)
    z = range(min_v, step = h, length = round(Int, (max_v - min_v) / h) + 1)

    return GridMatrix(x, y, z)
end