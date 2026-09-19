using QuantumMechanics

using Symbolics 
using SparseArrays
using KrylovKit


function eigen_solver(op::AbstractArray{T}, ::Type{T}; howmany::Int = min(10, size(op, 1)), which::Symbol = :SR) where T <: Real
    dim = size(op, 1)
    dim <= 0 && throw(DimensionMismatch("Hamiltonian dimension must be positive."))

    max_howmany = min(30, dim)
    if howmany > max_howmany
        @warn "Requested howmany=$howmany exceeds the stable KrylovKit limit for this matrix; reducing to $max_howmany."
        howmany = max_howmany
    end

    howmany = min(max(1, howmany), dim)
    v = Vector{T}(undef, dim)
    fill!(v, zero(T))
    v[1] = one(T)

    sol = eigsolve(op, v, howmany, which)
    eigenvals = T.(real.(sol[1]))
    eigen_vects = [T.(real.(e)) for e in sol[2]]

    return eigenvals, eigen_vects
end
