using QuantumMechanics
using LinearAlgebra
# Hydrogen ground-state example using the discretized Hamiltonian on a GridMatrix.
# This is a numerical approximation of the hydrogen atom problem, not the exact
# analytic solution.

T = Float64
L = T(12.0)
h = T(0.2)
x = range((-L + h/2), (L - h/2), step=h)
@time grid =GridMatrix(x, x, x)

# Hydrogen nucleus at the origin, with one electron charge.
@time V = CoulombPotential(-one(T), [one(T)], zeros(T, 1, 3))
@time H = hamiltonian(V, grid, Val(:central), T)

@time vals, vecs = eigen_solver(H, T; howmany = 2, which = :SR)
vals_sorted = sort(real.(vals))

for i in eachindex(vals)
    λ = vals[i]
    ψ = vecs[i]

    residual = norm(H * ψ - λ * ψ)

    println("λ[$i] = $λ")
    println("residual[$i] = $residual")
end

println("Hydrogen discretized spectrum:")
@show vals_sorted
println("Approximate ground-state energy:")
@show vals_sorted[1]
println("Corresponding eigenvector length:")
@show length(vecs[1])




