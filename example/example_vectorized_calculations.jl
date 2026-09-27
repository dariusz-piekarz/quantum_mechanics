using QuantumMechanics
using Symbolics

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