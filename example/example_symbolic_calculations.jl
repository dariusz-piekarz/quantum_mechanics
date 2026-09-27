using QuantumMechanics
using Symbolics


# ============================================================
# Symbolic calculations
# ============================================================

n, l, m = 2, 1, 1

@variables x α r θ φ

println("\n--- Symbolic ---")

@time P = legendre(l, x)
println("Legendre polynomial:")
display(P)

@time P_lm = associated_legendre(l, abs(m), x)
println("\nAssociated Legendre polynomial:")
display(P_lm)

@time L = generalized_laguerre(n - l - 1, 2l + 1, α)
println("\nGeneralized Laguerre polynomial:")
display(L)

@time Y = spherical_harmonic(l, m, θ, φ)
println("\nSpherical harmonic:")
display(Y)

@time ψ = wave_function(n, l, m, r, θ, φ)
println("\nHydrogen wave function:")
display(ψ)
