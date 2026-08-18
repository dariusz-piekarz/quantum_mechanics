using Symbolics
using QuantumMechanics

include("bench_benchmark.jl")

# ============================================================
# Symbolic variables
# ============================================================

@variables x y z

variables = (x, y, z)

expr = x^2 + y^2 + z^2

# ============================================================
# SymbolicFunction
# ============================================================

f = SymbolicFunction(expr, variables)

# ============================================================
# Precomputed symbolic expressions
# ============================================================

gradient_expr = gradient(expr, variables)
laplacian_expr = laplacian(expr, variables)
hessian_expr = hessian(expr, variables)

# ============================================================
# Precomputed SymbolicFunctions
# ============================================================

gradient_f = gradient(f)
laplacian_f = laplacian(f)
hessian_f = hessian(f)

# ============================================================
# DifferentialOperator
# ============================================================

run_benchmark(
    "D(x)(expr)",
    () -> D(x)(expr)
)

# ============================================================
# Gradient
# ============================================================

run_benchmark(
    "gradient(expr, variables)",
    () -> gradient(expr, variables)
)

run_benchmark(
    "gradient(f)",
    () -> gradient(f)
)

# ============================================================
# Laplacian
# ============================================================

run_benchmark(
    "laplacian(expr, variables)",
    () -> laplacian(expr, variables)
)

run_benchmark(
    "laplacian(f)",
    () -> laplacian(f)
)

# ============================================================
# Hessian
# ============================================================

run_benchmark(
    "hessian(expr, variables)",
    () -> hessian(expr, variables)
)

run_benchmark(
    "hessian(f)",
    () -> hessian(f)
)

# ============================================================
# Numerical evaluation
# ============================================================

run_benchmark(
    "f(1.0, 2.0, 3.0)",
    () -> f(1.0, 2.0, 3.0)
)

run_benchmark(
    "gradient_f(1.0, 2.0, 3.0)",
    () -> gradient_f(1.0, 2.0, 3.0)
)

run_benchmark(
    "laplacian_f(1.0, 2.0, 3.0)",
    () -> laplacian_f(1.0, 2.0, 3.0)
)

run_benchmark(
    "hessian_f(1.0, 2.0, 3.0)",
    () -> hessian_f(1.0, 2.0, 3.0)
)

println("\n", "=" ^ 70)
println("Benchmarking finished.")
println("=" ^ 70)