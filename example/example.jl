using QuantumMechanics
using Symbolics

# ============================================================
# Hydrogen atom
# ============================================================

n, l, m = 2, 1, 1
probability_level = Float32(0.1)

println("Hydrogen atom")
println("Quantum numbers: n = $n, l = $l, m = $m")


# ============================================================
# Symbolic calculations
# ============================================================

@variables x α r θ φ

println("\n--- Symbolic ---")

P = legendre(l, x)
println("Legendre polynomial:")
display(P)

P_lm = associated_legendre(l, abs(m), x)
println("\nAssociated Legendre polynomial:")
display(P_lm)

L = generalized_laguerre(n - l - 1, 2l + 1, α)
println("\nGeneralized Laguerre polynomial:")
display(L)

Y = spherical_harmonic(l, m, θ, φ)
println("\nSpherical harmonic:")
display(Y)

ψ = wave_function(n, l, m, r, θ, φ)
println("\nHydrogen wave function:")
display(ψ)


# ============================================================
# Numerical calculations
# ============================================================

println("\n--- Numerical ---")

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


# ============================================================
# Vectorized calculations
# ============================================================

println("\n--- Vectorized ---")

x_values = [0.5, 0.55]

P_values = associated_legendre(
    l,
    abs(m),
    x_values
)

println("Associated Legendre:")
display(P_values)

θ_values = [0.2, 0.4, 0.8]
φ_values = [3.1, 1.0, 2.0]

Y_values = spherical_harmonic(
    l,
    m,
    θ_values,
    φ_values
)

println("\nSpherical harmonics:")
display(Y_values)


# ============================================================
# Visualization
# ============================================================

println("\n--- Visualization ---")

# Legendre polynomial
display(plot_legendre(5))

# Associated Legendre polynomial
display(plot_associated_legendre(4, 2))

# Laguerre polynomial
display(plot_laguerre(9))

# Generalized Laguerre polynomial
display(plot_generalized_laguerre(9, 1.5))

# Spherical harmonic
display(
    plot_spherical_harmonic(
        l,
        m;
        draw=:real
    )
)

# Probability density
display(
    plot_probability_density(
        n,
        l,
        m;
        level=probability_level
    )
)

# Orbital
display(
    plot_orbital(
        n,
        l,
        m;
        level=probability_level
    )
)