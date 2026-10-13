using QuantumMechanics
using GLMakie, Meshing, GeometryBasics, LaTeXStrings

l, m = 3, 1

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
