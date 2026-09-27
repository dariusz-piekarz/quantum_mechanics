using QuantumMechanics
using Symbolics


# ============================================================
# Numerical calculations
# ============================================================

println("\n--- Numerical ---")

n, l, m = 2, 1, 1
x₀ = 0.2
α₀ = 0.4

println("P_$l($x₀) = ", legendre(l, x₀))

println(
    "P_$l^$(abs(m))($x₀) = ",
    associated_legendre(l, abs(m), x₀)
)

println(
    "L_$(n-l-1)^$(2l+1)($α₀) = ",
    generalized_laguerre(n - l - 1, 2l + 1, α₀)
)

println(
    "Y_$l^$m(0.2, 3.1) = ",
    spherical_harmonic(l, m, 0.2, 3.1)
)

println(
    "ψ($n, $l, $m; 0.1, 0.34, 0.66) = ",
    wave_function(n, l, m, 0.1, 0.34, 0.66)
)
