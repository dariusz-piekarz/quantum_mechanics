using BenchmarkTools
using Symbolics
using QuantumMechanics
using Printf

# ============================================================
# Symbolic function
# ============================================================

@variables x y z

expr = x^2 + sin(y) + exp(z)

f = SymbolicFunction(
    expr,
    (x, y, z)
)

# ============================================================
# Scalar evaluation
# ============================================================

function evaluate_scalar(f, X)

    npoints = size(X, 1)

    out = Vector{Float64}(undef, npoints)

    @inbounds for i in axes(X, 1)

        out[i] = f(
            X[i, 1],
            X[i, 2],
            X[i, 3]
        )

    end

    return out
end


# ============================================================
# Scalar evaluation with preallocated output
# ============================================================

function evaluate_scalar!(f, X, out)

    @inbounds for i in axes(X, 1)

        out[i] = f(
            X[i, 1],
            X[i, 2],
            X[i, 3]
        )

    end

    return out
end


# ============================================================
# Benchmark
# ============================================================

function run_benchmark()

    POINT_COUNTS = (
        1,
        10,
        100,
        1_000,
        10_000,
        100_000
    )

    points = Dict(
        n => rand(Float64, n, 3)
        for n in POINT_COUNTS
    )

    println("="^80)
    println("SCALAR vs BATCH EVALUATION")
    println("="^80)

    # ========================================================
    # Correctness
    # ========================================================

    println()
    println("CORRECTNESS")
    println("-"^80)

    X = points[100]

    scalar_result = evaluate_scalar(
        f.function_object,
        X
    )

    scalar_out = Vector{Float64}(undef, 100)

    scalar_result! = evaluate_scalar!(
        f.function_object,
        X,
        scalar_out
    )

    batch_result = evaluate(
        f,
        X
    )

    batch_out = Vector{Float64}(undef, 100)

    batch_result! = evaluate!(
        f,
        X,
        batch_out
    )

    @assert scalar_result ≈ batch_result
    @assert scalar_result ≈ scalar_result!
    @assert scalar_result ≈ batch_result!

    @assert scalar_result! === scalar_out
    @assert batch_result! === batch_out

    println("Correctness: OK")

    # ========================================================
    # Benchmarks
    # ========================================================

    for n in POINT_COUNTS

        X = points[n]

        scalar_out = Vector{Float64}(undef, n)
        batch_out  = Vector{Float64}(undef, n)

        println()
        println("-"^80)
        println("Points: $n")
        println("-"^80)

        # ----------------------------------------------------
        # Scalar
        # ----------------------------------------------------

        scalar_trial = @benchmark evaluate_scalar(
            $f.function_object,
            $X
        )

        # ----------------------------------------------------
        # Scalar!
        # ----------------------------------------------------

        scalar_preallocated_trial = @benchmark evaluate_scalar!(
            $f.function_object,
            $X,
            $scalar_out
        )

        # ----------------------------------------------------
        # Batch
        # ----------------------------------------------------

        batch_trial = @benchmark evaluate(
            $f,
            $X
        )

        # ----------------------------------------------------
        # Batch!
        # ----------------------------------------------------

        batch_preallocated_trial = @benchmark evaluate!(
            $f,
            $X,
            $batch_out
        )

        # ----------------------------------------------------
        # Median
        # ----------------------------------------------------

        scalar_time =
            median(scalar_trial).time

        scalar_preallocated_time =
            median(scalar_preallocated_trial).time

        batch_time =
            median(batch_trial).time

        batch_preallocated_time =
            median(batch_preallocated_trial).time

        # ----------------------------------------------------
        # Timing
        # ----------------------------------------------------

        println()

        @printf(
            "Scalar   : %s\n",
            BenchmarkTools.prettytime(scalar_time)
        )

        @printf(
            "Scalar!  : %s\n",
            BenchmarkTools.prettytime(
                scalar_preallocated_time
            )
        )

        @printf(
            "Batch    : %s\n",
            BenchmarkTools.prettytime(batch_time)
        )

        @printf(
            "Batch!   : %s\n",
            BenchmarkTools.prettytime(
                batch_preallocated_time
            )
        )

        # ----------------------------------------------------
        # Speedup
        # ----------------------------------------------------

        println()

        @printf(
            "Batch / Scalar       : %.2fx\n",
            batch_time / scalar_time
        )

        @printf(
            "Batch! / Scalar!     : %.2fx\n",
            batch_preallocated_time /
            scalar_preallocated_time
        )

        @printf(
            "Batch / Scalar!      : %.2fx\n",
            batch_time /
            scalar_preallocated_time
        )

        @printf(
            "Batch! / Scalar      : %.2fx\n",
            batch_preallocated_time /
            scalar_time
        )

        # ----------------------------------------------------
        # Allocations
        # ----------------------------------------------------

        println()

        @printf(
            "Scalar allocs        : %d\n",
            scalar_trial.allocs
        )

        @printf(
            "Scalar! allocs       : %d\n",
            scalar_preallocated_trial.allocs
        )

        @printf(
            "Batch allocs         : %d\n",
            batch_trial.allocs
        )

        @printf(
            "Batch! allocs        : %d\n",
            batch_preallocated_trial.allocs
        )

        # ----------------------------------------------------
        # Memory
        # ----------------------------------------------------

        println()

        @printf(
            "Scalar memory        : %s\n",
            Base.format_bytes(
                scalar_trial.memory
            )
        )

        @printf(
            "Scalar! memory       : %s\n",
            Base.format_bytes(
                scalar_preallocated_trial.memory
            )
        )

        @printf(
            "Batch memory         : %s\n",
            Base.format_bytes(
                batch_trial.memory
            )
        )

        @printf(
            "Batch! memory        : %s\n",
            Base.format_bytes(
                batch_preallocated_trial.memory
            )
        )
    end

    println()
    println("="^80)
    println("BENCHMARK FINISHED")
    println("="^80)

end


run_benchmark()