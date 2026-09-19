using QuantumMechanics

# Hydrogen ground-state example using the discretized Hamiltonian on a GridMatrix.
# This is a numerical approximation of the hydrogen atom problem, not the exact
# analytic solution.

T = Float64
x = range(0.1, step = 0.1, length = 10)
grid = GridMatrix(x, x, x)

# Hydrogen nucleus at the origin, with one electron charge.
V = CoulombPotential(-one(T), [one(T)], zeros(T, 1, 3))
H = hamiltonian(V, grid, Val(:central), T)

vals, vecs = eigen_solver(H, T; howmany = 3, which = :SR)
vals_sorted = sort(real.(vals))

println("Hydrogen discretized spectrum:")
@show vals_sorted
println("Approximate ground-state energy:")
@show vals_sorted[1]
println("Corresponding eigenvector length:")
@show length(vecs[1])


