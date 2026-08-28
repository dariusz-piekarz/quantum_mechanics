using BenchmarkTools
using Printf


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

