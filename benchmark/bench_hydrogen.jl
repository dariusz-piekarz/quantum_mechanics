using QuantumMechanics


include("bench_benchmark.jl")

# ============================================================
# Data
# ============================================================


const ξ20 = collect(range(-10.0, 10.0, length=20))
const ξ50 = collect(range(-10.0, 10.0, length=50))

r = rand(Float32, 10_000)
const θ_values = rand(Float64, 10_000) .* π
const φ_values = (rand(Float64, 10_000) .- 0.5) .* 2π

# ============================================================
# Hydrogen wave function
# ============================================================

wave_trial = run_benchmark(
"Hydrogen wave function",
    () -> wave_function(3, 2, 1, r, θ_values, φ_values)
)


# ============================================================
# Cartesian grid 20³
# ============================================================

grid20_trial = run_benchmark(
    "Cartesian wave function 20³",
    () -> wave_function_cartesian_grid(3, 2, 1, ξ20)
)


# ============================================================
# Cartesian grid 50³
# ============================================================

grid50_trial = run_benchmark(
    "Cartesian wave function 50³",
    () -> wave_function_cartesian_grid(3, 2, 1, ξ50)
)


println("\n", "=" ^ 70)
println("Benchmarking finished.")
println("=" ^ 70)