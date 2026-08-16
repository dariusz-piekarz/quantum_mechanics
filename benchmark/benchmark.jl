using BenchmarkTools
using QuantumMechanics
using Printf

# ============================================================
# Data
# ============================================================

x = collect(range(-1.0, 1.0, length=10_000))

θ = collect(range(0.0, π, length=100))
φ = collect(range(0.0, 2π, length=100))

const ξ20 = collect(range(-10.0, 10.0, length=20))
const ξ50 = collect(range(-10.0, 10.0, length=50))

r = rand(Float32, 10_000)
const θ_values = rand(Float64, 10_000) .* π
const φ_values = (rand(Float64, 10_000) .- 0.5) .* 2π


# ============================================================
# Benchmark helper
# ============================================================

function run_benchmark(name::String, f::Function)
    println("\n", "-" ^ 70)
    println(name)
    println("-" ^ 70)

    trial = @benchmark $f()

    display(trial)

    println()
    println("Median time : ", median(trial))
    println("Minimum time: ", minimum(trial))
    println("Maximum time: ", maximum(trial))
    println("Allocations : ", trial.allocs)
    println("Memory      : ", Base.format_bytes(trial.memory))

    return trial
end


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