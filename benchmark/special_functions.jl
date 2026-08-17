using QuantumMechanics

include("benchmark.jl")

# ============================================================
# Data
# ============================================================

x = collect(range(-1.0, 1.0, length=10_000))

θ = collect(range(0.0, π, length=100))
φ = collect(range(0.0, 2π, length=100))

# ============================================================
# Associated Legendre
# ============================================================

legendre_trial = run_benchmark(
    "Associated Legendre P₅²",
    () -> associated_legendre(5, 2, x)
)


# ============================================================
# Generalized Laguerre
# ============================================================

laguerre_trial = run_benchmark(
    "Generalized Laguerre L₅²",
    () -> generalized_laguerre(5, 2, x)
)


# ============================================================
# Spherical harmonics
# ============================================================

spherical_trial = run_benchmark(
    "Spherical harmonic Y₅²",
    () -> spherical_harmonic(5, 2, θ, φ)
)


println("\n", "=" ^ 70)
println("Benchmarking finished.")
println("=" ^ 70)
