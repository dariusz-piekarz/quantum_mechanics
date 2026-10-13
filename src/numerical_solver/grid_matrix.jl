using QuantumMechanics


"""
    GridMatrix(x, y, z)
    GridMatrix(min, max, h)
    GridMatrix(xmin, xmax, ymin, ymax, zmin, zmax, hx, hy, hz)

Cartesian grid with `x` varying fastest, followed by `y` and `z`. Each axis
needs at least two points and may have its own spacing. For bounds and a step,
`min` is included and `max` is an inclusive upper limit: the last point is
`min + k*h ≤ max`, so `max` need not be a grid point.

Coordinates and steps have the same units. The numerical Hamiltonian uses
atomic units: provide coordinates in bohr and potentials in hartree to obtain
energies in hartree. No unit conversion is performed.

The finite-difference operator uses zero-valued ghost points one step beyond
the first and last sample on each axis (homogeneous Dirichlet boundaries).
"""
struct GridMatrix{T<:Real}
    grid::Matrix{T}

    function GridMatrix{T}(grid::AbstractMatrix{<:Real}) where {T <: Real}
        size(grid, 2) == 3 || throw(DimensionMismatch("GridMatrix coordinates must have three columns (x, y, z)."))
        return new{T}(Matrix{T}(grid))
    end
end


function GridMatrix(x::AbstractRange{<:Real}, y::AbstractRange{<:Real}, z::AbstractRange{<:Real})
    if length(x) < 2 || length(y) < 2 || length(z) < 2
        throw(DimensionMismatch("Grid dimensions must contain at least two points."))
    end

    for axis in (x, y, z)
        (isfinite(first(axis)) && isfinite(last(axis)) && isfinite(step(axis)) && !iszero(step(axis))) ||
            throw(ArgumentError("Grid coordinates and steps must be finite, with nonzero steps."))
    end

    T = promote_type(typeof(float(first(x))), typeof(float(first(y))), typeof(float(first(z))))
    nx, ny, nz = length(x), length(y), length(z)
    grid = Matrix{T}(undef, nx * ny * nz, 3)
    idx = 1
    @inbounds for k in z, j in y, i in x
        grid[idx, 1] = i
        grid[idx, 2] = j
        grid[idx, 3] = k
        idx += 1
    end

    return GridMatrix{T}(grid)
end


GridMatrix(base_range::AbstractRange{<:Real}) = GridMatrix(base_range, base_range, base_range)


function _bounded_axis(min_v::Real, max_v::Real, h::Real)
    min_v, max_v, h = promote(float(min_v), float(max_v), float(h))
    all(isfinite, (min_v, max_v, h)) ||
        throw(ArgumentError("Grid bounds and step must be finite."))
    max_v > min_v || throw(DimensionMismatch("Grid max must be greater than min."))
    h > zero(h) || throw(DimensionMismatch("Grid step must be positive."))

    # A stop-based range never rounds the number of steps up past max_v.
    axis = range(min_v; step=h, stop=max_v)
    length(axis) >= 2 ||
        throw(DimensionMismatch("Grid bounds and step must yield at least two points per axis."))
    return axis
end


function GridMatrix(
    x_min::Real, x_max::Real,
    y_min::Real, y_max::Real,
    z_min::Real, z_max::Real,
    hx::Real, hy::Real, hz::Real,
)
    return GridMatrix(
        _bounded_axis(x_min, x_max, hx),
        _bounded_axis(y_min, y_max, hy),
        _bounded_axis(z_min, z_max, hz),
    )
end


GridMatrix(x_min::Real, x_max::Real, y_min::Real, y_max::Real,
           z_min::Real, z_max::Real, h::Real) =
    GridMatrix(x_min, x_max, y_min, y_max, z_min, z_max, h, h, h)


GridMatrix(min_v::Real, max_v::Real, h::Real) =
    GridMatrix(min_v, max_v, min_v, max_v, min_v, max_v, h)
