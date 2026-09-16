using BenchmarkTools
using Printf
using Symbolics
using QuantumMechanics

# Compare direct Cartesian iteration over Grid with batch evaluation over
# the point matrix owned by GridMatrix.
const T = Float64
const N = 32
const STEP = T(0.25)

@variables x y z

const ψ = SymbolicFunction(x^2 + y^2 + z^2, (x, y, z))
const H_sym = Hamiltonian(HarmonicPotential(one(T)), KineticOperator(T), ψ)
const axis = range(-T(4), step = STEP, length = N)
const grid = Grid(axis, axis, axis)
const grid_matrix = GridMatrix(axis, axis, axis)

function evaluate_grid!(f::SymbolicFunction, grid::Grid, output::AbstractVector)
    point = 1
    @inbounds for x in grid.x, y in grid.y, z in grid.z
        output[point] = f.function_object(x, y, z)
        point += 1
    end
    return output
end

function evaluate_grid(f::SymbolicFunction, grid::Grid)
    output = Vector{T}(undef, length(grid.x) * length(grid.y) * length(grid.z))
    return evaluate_grid!(f, grid, output)
end

function evaluate_grid_matrix(f::SymbolicFunction, grid::GridMatrix)
    return evaluate(f, grid.grid)
end

grid_output = Vector{T}(undef, N^3)
matrix_output = Vector{T}(undef, N^3)

evaluate_grid!(H_sym, grid, grid_output)
evaluate!(H_sym, grid_matrix.grid, matrix_output)
@assert grid_output ≈ matrix_output

function report(name::String, trial)
    estimate = median(trial)
    @printf(
        "%-28s %12s %12s %8d\n",
        name,
        BenchmarkTools.prettytime(estimate.time),
        Base.format_bytes(estimate.memory),
        estimate.allocs
    )
end

println("Hamiltonian evaluation: Grid vs GridMatrix")
println("Grid points: ", N, "^3 = ", N^3)
println("Symbolic expression: ", H_sym.expression)
println()
@printf("%-28s %12s %12s %8s\n", "Case", "Median", "Memory", "Allocs")
println("-"^68)

grid_trial = @benchmark evaluate_grid!($H_sym, $grid, $grid_output)
matrix_trial = @benchmark evaluate!($H_sym, $grid_matrix.grid, $matrix_output)
grid_allocating_trial = @benchmark evaluate_grid($H_sym, $grid)
matrix_allocating_trial = @benchmark evaluate_grid_matrix($H_sym, $grid_matrix)
grid_matrix_creation_trial = @benchmark GridMatrix($grid.x, $grid.y, $grid.z)

report("Grid, preallocated", grid_trial)
report("GridMatrix, preallocated", matrix_trial)
report("Grid, allocating", grid_allocating_trial)
report("GridMatrix, allocating", matrix_allocating_trial)
report("GridMatrix construction", grid_matrix_creation_trial)

println()
@printf(
    "Preallocated speedup (Grid/Matrix): %.2fx\n",
    median(grid_trial).time / median(matrix_trial).time
)
@printf(
    "Preallocated memory ratio (Grid/Matrix): %.2fx\n",
    median(grid_trial).memory / max(median(matrix_trial).memory, 1)
)