using QuantumMechanics

using SparseArrays
using KrylovKit
using LinearAlgebra: issymmetric, eigen, Symmetric


function eigen_solver(op::AbstractMatrix{T}, ::Type{T}; howmany::Int = min(10, size(op, 1)), which::Symbol = :SR) where T <: AbstractFloat
    dim = size(op, 1)
    dim <= 0 && throw(DimensionMismatch("Hamiltonian dimension must be positive."))
    size(op, 2) == dim || throw(DimensionMismatch("Hamiltonian must be square."))
    1 <= howmany <= dim || throw(ArgumentError("howmany must be between 1 and the Hamiltonian dimension."))
    issymmetric(op) || throw(ArgumentError("The Hamiltonian must be symmetric for the Lanczos eigensolver."))

    # For small matrices, a dense symmetric solve is more reliable than block
    # Lanczos near an exhausted Krylov subspace, particularly at degeneracies.
    if dim <= 128 && which in (:SR, :LR, :LM)
        decomposition = eigen(Symmetric(Matrix(op)))
        order = if which == :SR
            sortperm(decomposition.values)
        elseif which == :LR
            sortperm(decomposition.values; rev=true)
        else
            sortperm(decomposition.values; by=abs, rev=true)
        end
        indices = order[1:howmany]
        return decomposition.values[indices], [decomposition.vectors[:, i] for i in indices]
    end

    # A single basis vector can be orthogonal to the desired eigenspace (for
    # example, e₁ in a diagonal matrix). A block also resolves degeneracies.
    initial_vectors = Block([rand(T, dim) for _ in 1:howmany])

    vals, vecs, info = eigsolve(op, initial_vectors, howmany, which, tol = T(1e-4), maxiter=150, krylovdim=max(30, 2howmany + 1), ishermitian=true, issymmetric=true, eager=true)
    info.converged >= howmany ||
        error("Eigensolver converged for $(info.converged) of $howmany requested eigenpairs.")
    # KrylovKit may return extra converged eigenpairs from the final block.
    eigenvals = T.(real.(vals[1:howmany]))
    eigen_vects = [T.(real.(e)) for e in vecs[1:howmany]]

    return eigenvals, eigen_vects
end
