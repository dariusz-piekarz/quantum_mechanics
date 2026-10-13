using QuantumMechanics

using SparseArrays


function kinetic_operator(
    n::Integer,
    h::T,
    ::Type{T})::SparseMatrixCSC{T} where T <: Real
    return laplacian(n, h, T, 3, -one(T) / T(2))
end


function kinetic_operator(
    grid::GridMatrix{T},
    ::Type{T})::SparseMatrixCSC{T} where T <: Real
    n = round(Int, cbrt(size(grid.grid, 1)))

    n >= 2 && n^3 == size(grid.grid, 1) ||
        throw(DimensionMismatch("GridMatrix must represent a cubic grid with at least two points per axis."))

    coords = grid.grid
    x0, y0, z0 = coords[1, 1], coords[1, 2], coords[1, 3]
    h = coords[2, 1] - x0
    (iszero(h) || !isfinite(h)) && throw(ArgumentError("Grid spacing must be finite and nonzero."))

    # The Kronecker stencil assumes x varies fastest, then y, then z, with
    # equal spacing along all three axes. Checking only the total point count
    # would silently accept e.g. a 2 × 4 × 8 grid as a 4 × 4 × 4 grid.
    tol = 8n * maximum(eps(float(v)) for v in (x0, y0, z0, h))
    @inbounds for k in 1:n, j in 1:n, i in 1:n
        idx = i + n * (j - 1 + n * (k - 1))
        if !(isapprox(coords[idx, 1] - x0, (i - 1) * h; atol=tol, rtol=0) &&
             isapprox(coords[idx, 2] - y0, (j - 1) * h; atol=tol, rtol=0) &&
             isapprox(coords[idx, 3] - z0, (k - 1) * h; atol=tol, rtol=0))
            throw(DimensionMismatch("GridMatrix must be a regular cubic grid with equal spacing in x, y, z and x varying fastest."))
        end
    end

    return laplacian(n, h, T, 3, -one(T) / T(2))
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
