using Test
using Symbolics
using SparseArrays
using QuantumMechanics

@testset "Grid" begin
    g = Grid(3, 0.1, Float64)

    @test g isa Grid{Float64}
    @test grid_size(g) == (3, 3, 3)
    @test grid_step(g) == 0.1
    @test length(g.x) == 3
    @test length(g.y) == 3
    @test length(g.z) == 3
end

@testset "Grid validation" begin
    @test_throws DimensionMismatch Grid(0, 0.1, Float64)

    x = range(0.0, step = 0.1, length = 3)
    y = range(0.0, step = 0.2, length = 3)
    z = range(0.0, step = 0.1, length = 3)

    @test_throws DimensionMismatch Grid(x, y, z)
end

@testset "GridMatrix" begin
    T = Float64
    x = range(0.0, step = 0.1, length = 3)
    y = range(0.0, step = 0.1, length = 3)
    z = range(0.0, step = 0.1, length = 3)

    gm = GridMatrix(x, y, z)

    @test gm isa GridMatrix{Float64}
    @test size(gm.grid) == (27, 3)
    @test gm.grid[1, :] == [0.0, 0.0, 0.0]
    @test gm.grid[end, :] == [0.2, 0.2, 0.2]

    @test_throws DimensionMismatch GridMatrix(range(0.0, step = 0.1, length = 1), y, z)
    @test_throws DimensionMismatch GridMatrix(x, range(0.0, step = 0.2, length = 3), z)
end

@testset "potential_operator" begin
    T = Float64
    g = Grid(2, 0.5, T)
    V = HarmonicPotential(1.0)
    theoretical_potential = vec([V(x,y,z) for x in g.x, y in g.y, z in g.z])
    R = spdiagm(0 => theoretical_potential)

    P = potential_operator(V)(g, T)

    @test P isa SparseMatrixCSC{T}
    @test size(P) == (8, 8)
    @test P == R

end

@testset "kinetic_operator" begin
    T = Float64
    g = Grid(2, 0.5, T)
    K = kinetic_operator(g, Val(:central), T)

    @test K isa SparseMatrixCSC{T}
    @test size(K) == (8, 8)
end

@testset "GridMatrix kinetic_operator and hamiltonian" begin
    T = Float64
    x = range(0.0, step = 0.5, length = 2)
    gm = GridMatrix(x, x, x)
    V = HarmonicPotential(1.0)

    K = kinetic_operator(gm, Val(:central), T)
    P = potential_operator(V)(gm, T)
    H = hamiltonian(V, gm, Val(:central), T)

    @test K isa SparseMatrixCSC{T}
    @test P isa SparseMatrixCSC{T}
    @test H isa SparseMatrixCSC{T}
    @test size(K) == (8, 8)
    @test size(P) == (8, 8)
    @test size(H) == (8, 8)
    @test H[1, 1] ≈ 12.0
    @test H[2, 2] ≈ 12.125
    @test H[8, 8] ≈ 12.375
end

@testset "hamiltonian" begin
    T = Float64
    g = Grid(2, 0.5, T)
    V = HarmonicPotential(1.0)

    H = hamiltonian(V, g, Val(:central), T)

    @test H isa SparseMatrixCSC{T}
    @test size(H) == (8, 8)

    # check diagonal consistency for the harmonic potential
    @test H[1, 1] ≈ 12.0
    @test H[2, 2] ≈ 12.125
    @test H[8, 8] ≈ 12.375
end


@testset "symbolic Hamiltonian matches numerical Hamiltonian on the interior grid" begin
    T = Float64
    @variables x y z

    ψ_expr = x^2 + y^2 + z^2
    ψ = SymbolicFunction(ψ_expr, (x, y, z))
    K = KineticOperator(T)

    x_axis = range(0.5, step = 0.5, length = 3)
    grid = Grid(x_axis, x_axis, x_axis)

    potentials = [
        ("V_e", V_e(T)),
        ("V_ee", V_ee(T)),
        ("V_en", V_en(T)),
        ("V_ne", V_ne(T)),
        ("Harmonic", HarmonicPotential(T(1.0)))
    ]

    for (name, V) in potentials
        H_sym = Hamiltonian(V, K, ψ)
        H_num = hamiltonian(V, grid, Val(:central), T)

        ψ_values = T[
            ψ.function_object(x, y, z)
            for x in grid.x, y in grid.y, z in grid.z
        ]

        H_sym_values = T[
            H_sym.function_object(x, y, z)
            for x in grid.x, y in grid.y, z in grid.z
        ]

        # The discrete finite-difference Laplacian is only expected to match the
        # symbolic operator on the interior of the grid. The boundary stencil is
        # a truncated approximation and therefore differs by construction.
        H_num_values = reshape(H_num * vec(ψ_values), size(ψ_values))

        H_num_interior = H_num_values[2:end-1, 2:end-1, 2:end-1]
        H_sym_interior = H_sym_values[2:end-1, 2:end-1, 2:end-1]

        @test H_num_interior ≈ H_sym_interior atol = 1e-10 rtol = 1e-10
    end
end


@testset "symbolic Hamiltonian matches numerical Hamiltonian on the interior grid (matrix based)" begin
    T = Float64
    @variables x y z

    ψ_expr = x^2 + y^2 + z^2
    ψ = SymbolicFunction(ψ_expr, (x, y, z))
    K = KineticOperator(T)

    x_axis = range(0.5, step = 0.5, length = 3)
    grid = GridMatrix(x_axis, x_axis, x_axis)
    ψ_values = ψ(grid.grid)

    potentials = [
        ("V_e", V_e(T)),
        ("V_ee", V_ee(T)),
        ("V_en", V_en(T)),
        ("V_ne", V_ne(T)),
        ("Harmonic", HarmonicPotential(T(1.0)))
    ]

    for (name, V) in potentials
        H_sym = Hamiltonian(V, K, ψ)
        H_num = hamiltonian(V, grid, Val(:central), T)

        H_sym_values = Vector{T}(undef, size(grid.grid, 1))
        evaluate!(ψ, grid.grid, H_sym_values)

        # The discrete finite-difference Laplacian is only expected to match the
        # symbolic operator on the interior of the grid. The boundary stencil is
        # a truncated approximation and therefore differs by construction.
        H_num_values = reshape(H_num * vec(ψ_values), size(ψ_values))

        H_num_interior = H_num_values[2:end-1, 2:end-1, 2:end-1]
        H_sym_interior = H_sym_values[2:end-1, 2:end-1, 2:end-1]

        @test H_num_interior ≈ H_sym_interior atol = 1e-10 rtol = 1e-10
    end
end
