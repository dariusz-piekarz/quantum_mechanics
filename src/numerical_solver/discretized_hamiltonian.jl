using QuantumMechanics

using SparseArrays


function kinetic_operator(
    n::Integer,
    h::T,
    ::Type{T})::SparseMatrixCSC{T} where T <: Real
    return laplacian(n, h, T, 3, -one(T) / T(2))
end


function _grid_geometry(grid::GridMatrix{T}) where T <: Real
    coords = grid.grid
    npoints = size(coords, 1)
    npoints >= 8 || throw(DimensionMismatch("Grid needs at least two points per axis."))

    x0, y0, z0 = coords[1, 1], coords[1, 2], coords[1, 3]
    all(isfinite, (x0, y0, z0)) ||
        throw(ArgumentError("Grid coordinates must be finite."))

    # Find the first y and z transitions in the flattened x-fastest layout.
    nx = 1
    while nx < npoints && coords[nx + 1, 2] == y0 && coords[nx + 1, 3] == z0
        nx += 1
    end
    nxy = 1
    while nxy < npoints && coords[nxy + 1, 3] == z0
        nxy += 1
    end
    nx < npoints && nxy < npoints && nxy % nx == 0 && npoints % nxy == 0 ||
        throw(DimensionMismatch("Grid must be a regular Cartesian product with x varying fastest."))

    ny, nz = nxy ÷ nx, npoints ÷ nxy
    all(n -> n >= 2, (nx, ny, nz)) ||
        throw(DimensionMismatch("Grid needs at least two points per axis."))

    hx = coords[2, 1] - x0
    hy = coords[nx + 1, 2] - y0
    hz = coords[nxy + 1, 3] - z0
    steps = (hx, hy, hz)
    all(h -> isfinite(h) && !iszero(h), steps) ||
        throw(ArgumentError("Grid spacing must be finite and nonzero on every axis."))

    origins = (x0, y0, z0)
    dims = (nx, ny, nz)
    tolerances = ntuple(3) do axis
        origin, h, n = origins[axis], steps[axis], dims[axis]
        extent = origin + (n - 1) * h
        scale = max(abs(origin), abs(extent), abs(h))
        min(8n * eps(float(scale)), abs(h) / 4)
    end

    @inbounds for k in 1:nz, j in 1:ny, i in 1:nx
        idx = i + nx * (j - 1 + ny * (k - 1))
        if !(isapprox(coords[idx, 1], x0 + (i - 1) * hx; atol=tolerances[1], rtol=0) &&
             isapprox(coords[idx, 2], y0 + (j - 1) * hy; atol=tolerances[2], rtol=0) &&
             isapprox(coords[idx, 3], z0 + (k - 1) * hz; atol=tolerances[3], rtol=0))
            throw(DimensionMismatch("Grid must be regular, with x varying fastest."))
        end
    end

    return dims, steps
end


function kinetic_operator(
    grid::GridMatrix{T},
    ::Type{T})::SparseMatrixCSC{T} where T <: Real
    T <: AbstractFloat || throw(ArgumentError("Grid coordinates must use a floating-point type for the kinetic operator."))
    dims, steps = _grid_geometry(grid)

    return laplacian(dims, steps, T, -one(T) / T(2))
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


function (op::PotentialOperator)(grid::GridMatrix{T}, ::Type{T})::SparseMatrixCSC{T} where {T <: Real}
    diag = op.potential(grid)

    return spdiagm(0 => diag)
end


function potential_operator(V::AbstractPotential)
    return PotentialOperator(V)
end


function hamiltonian(
    V::AbstractPotential,
    grid::GridMatrix{T},
    ::Type{T}
)::SparseMatrixCSC{T} where T <: Real

    return kinetic_operator(grid, T) + potential_operator(V)(grid, T)
end
