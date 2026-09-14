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


function kinetic_operator(
    grid::GridMatrix{T},
    cut_type::Val,
    ::Type{T})::SparseMatrixCSC{T} where T <: Real
    n = round(Int, cbrt(size(grid.grid, 1)))

    if n^3 != size(grid.grid, 1)
        throw(DimensionMismatch("GridMatrix must represent a cubic grid for the 3D discretized Laplacian."))
    end

    return -one(T) / T(2) * laplacian(n, cut_type, grid.grid[2,1] - grid.grid[1,1], T, 3)
end


struct PotentialOperator{VType <: AbstractPotential}
    potential::VType
end

function (V::AbstractPotential)(grid::GridMatrix{T})::Vector{T} where {T <: Real}
    diag = Vector{T}(undef, size(grid.grid, 1))

    @inbounds for i in 1:size(grid.grid, 1)
        diag[i] = V(grid.grid[i, 1], grid.grid[i, 2], grid.grid[i, 3])
    end

    return diag
end

function (op::PotentialOperator)(grid::Grid{T}, ::Type{T})::SparseMatrixCSC{T} where {T <: Real}
    diag = Vector{T}(undef, length(grid.x) * length(grid.y) * length(grid.z))
    k = 1

    for x in grid.x, y in grid.y, z in grid.z
        diag[k] = op.potential(x, y, z)
        k += 1
    end

    return spdiagm(0 => diag)
end


function (op::PotentialOperator)(grid::GridMatrix{T}, ::Type{T})::SparseMatrixCSC{T} where {T <: Real}
    diag = op.potential(grid)

    return spdiagm(0 => diag)
end


function potential_operator(V::AbstractPotential)
    return PotentialOperator(V)
end


function hamiltonian(
    V::AbstractPotential,
    grid::Grid{T},
    cut_type::Union{Val{:forward}, Val{:backward}, Val{:central}},
    ::Type{T}
)::SparseMatrixCSC{T} where T <: Real

    return kinetic_operator(grid, cut_type, T) + potential_operator(V)(grid, T)
end


function hamiltonian(
    V::AbstractPotential,
    grid::GridMatrix{T},
    cut_type::Union{Val{:forward}, Val{:backward}, Val{:central}},
    ::Type{T}
)::SparseMatrixCSC{T} where T <: Real

    return kinetic_operator(grid, cut_type, T) + potential_operator(V)(grid, T)
end


function hamiltonian(
    V::AbstractPotential,
    n::Integer,
    cut_type::Union{Val{:forward}, Val{:backward}, Val{:central}},
    h::T,
    ::Type{T}
)::SparseMatrixCSC{T} where T <: Real

    return hamiltonian(V, Grid(n, h, T), cut_type, T)
end

