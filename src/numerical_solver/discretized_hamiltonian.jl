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

    if n^3 != size(grid.grid, 1)
        throw(DimensionMismatch("GridMatrix must represent a cubic grid for the 3D discretized Laplacian."))
    end

    return laplacian(n, grid.grid[2,1] - grid.grid[1,1], T, 3, -one(T) / T(2))
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
