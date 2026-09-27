using QuantumMechanics
using Symbolics

# ============================================================
# Hydrogen atom
# ============================================================

n, l, m = 2, 1, 1
probability_level = Float32(0.1)

println("\n--- Hydrogen atom ---")
println("Quantum numbers: n = $n, l = $l, m = $m")


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
