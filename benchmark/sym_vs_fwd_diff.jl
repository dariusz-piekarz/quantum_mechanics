using Symbolics
using ForwardDiff
using QuantumMechanics

include("benchmark.jl")

# ============================================================
# Configuration
# ============================================================

const N = 10

# ============================================================
# Symbolic variables
# ============================================================

@variables x1 x2 x3 x4 x5 x6 x7 x8 x9 x10

const variables = (
    x1, x2, x3, x4, x5,
    x6, x7, x8, x9, x10
)

# ============================================================
# Complex symbolic function
# ============================================================

expr =
    0.05 * sin(x1 * x2 * x3) +

    sum(
        exp(-0.5 * variables[i]^2) * sin(variables[i])
        for i in 1:N
    ) +

    sum(
        0.1 * variables[i] * variables[i + 1]
        for i in 1:N-1
    ) +

    sum(
        variables[i]^3 / 10
        for i in 1:N
    )

println("\nSymbolic expression:")
@show expr

# ============================================================
# Symbolics -> symbolic gradient
# ============================================================

println("\nGenerating symbolic gradient...")

@time gradient_expr = gradient(expr, variables)

# ============================================================
# Symbolics -> numerical function
# ============================================================

println("\nGenerating numerical function...")

@time gradient_function = to_function(
    gradient_expr,
    variables
)

# ============================================================
# Pure Julia function for ForwardDiff
# ============================================================

function f_forward(x)

    result =
        0.05 * sin(x[1] * x[2] * x[3])

    @inbounds for i in 1:10
        result +=
            exp(-0.5 * x[i]^2) *
            sin(x[i])
    end

    @inbounds for i in 1:9
        result +=
            0.1 * x[i] * x[i + 1]
    end

    @inbounds for i in 1:10
        result +=
            x[i]^3 / 10
    end

    return result
end

# ============================================================
# Test point
# ============================================================

const point = rand(Float64, N)

println("\nPoint:")
@show point

# ============================================================
# Correctness
# ============================================================

symbolic_result =
    gradient_function(point...)

forward_result =
    ForwardDiff.gradient(f_forward, point)

println("\nCorrectness check:")

@show symbolic_result
@show forward_result

@assert symbolic_result ≈ forward_result

println("\nGradient results agree.")

# ============================================================
# Benchmark: gradient construction
# ============================================================

println("\n")
println("=" ^ 70)
println("GRADIENT CONSTRUCTION")
println("=" ^ 70)

run_benchmark(
    "Symbolics: symbolic gradient construction",
    () -> gradient(expr, variables)
)

# ============================================================
# Benchmark: numerical gradient evaluation
# ============================================================

println("\n")
println("=" ^ 70)
println("NUMERICAL GRADIENT EVALUATION")
println("=" ^ 70)

run_benchmark(
    "Symbolics -> numerical function",
    () -> gradient_function(point...)
)

run_benchmark(
    "ForwardDiff",
    () -> ForwardDiff.gradient(f_forward, point)
)

# ============================================================
# Multiple points
# ============================================================

println("\n")
println("=" ^ 70)
println("MULTIPLE POINTS")
println("=" ^ 70)


function evaluate_symbolic!(
    output,
    f,
    points
)

    @inbounds for i in axes(points, 2)
        output[:, i] .= f(points[:, i]...)
    end

    return output
end


function evaluate_forward!(
    output,
    f,
    points
)

    @inbounds for i in axes(points, 2)

        output[:, i] .= ForwardDiff.gradient(
            f,
            view(points, :, i)
        )

    end

    return output
end


# ============================================================
# Different numbers of points
# ============================================================

for number_of_points in (
    1,
    100,
    1_000,
    10_000,
)

    println("\n")
    println("-" ^ 70)
    println("Points: ", number_of_points)
    println("-" ^ 70)

    points = rand(
        Float64,
        N,
        number_of_points
    )

    symbolic_output = Matrix{Float64}(
        undef,
        N,
        number_of_points
    )

    forward_output = Matrix{Float64}(
        undef,
        N,
        number_of_points
    )

    # --------------------------------------------------------
    # Correctness
    # --------------------------------------------------------

    evaluate_symbolic!(
        symbolic_output,
        gradient_function,
        points
    )

    evaluate_forward!(
        forward_output,
        f_forward,
        points
    )

    @assert symbolic_output ≈ forward_output

    # --------------------------------------------------------
    # Symbolics -> numerical function
    # --------------------------------------------------------

    run_benchmark(
        "Symbolics -> function ($number_of_points points)",
        () -> evaluate_symbolic!(
            symbolic_output,
            gradient_function,
            points
        )
    )

    # --------------------------------------------------------
    # ForwardDiff
    # --------------------------------------------------------

    run_benchmark(
        "ForwardDiff ($number_of_points points)",
        () -> evaluate_forward!(
            forward_output,
            f_forward,
            points
        )
    )
end

println("\n")
println("=" ^ 70)
println("Benchmarking finished.")
println("=" ^ 70)
