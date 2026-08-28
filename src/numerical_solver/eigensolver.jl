using QuantumMechanics

using Symbolics 
using SparseArrays
using KrylovKit


function eigen_solver(op::AbstractArray{T}, ::Type{T}) where T <: Real
    dim = size(op)[1]
    v = [(e > 1 ? zero(T) : one(T)) for e in 1: dim]
    sol = eigsolve(op, v, dim, :LM)
    eigenvals = T.(real.(sol[1]))
    eigen_vects = [T.(real.(e)) for e in sol[2]]

    return eigenvals, eigen_vects
end
