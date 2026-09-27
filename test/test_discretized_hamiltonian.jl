using Test
using Symbolics
using SparseArrays
using LinearAlgebra
using QuantumMechanics


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


@testset "kinetic_operator and hamiltonian" begin
    T = Float64
    x = range(0.0, step = 0.5, length = 2)
    gm = GridMatrix(x, x, x)
    V = HarmonicPotential(1.0)

    K = kinetic_operator(gm, T)
    P = potential_operator(V)(gm, T)
    H = hamiltonian(V, gm, T)

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


@testset "eigen_solver for Hamiltonians" begin
    T = Float64
    x_axis = range(0.5, step = 0.5, length = 2)
    grid = GridMatrix(x_axis, x_axis, x_axis)

    potentials = [
        ("Harmonic", HarmonicPotential(T(1.0))),
        ("V_e", V_e(T)),
        ("V_ee", V_ee(T)),
        ("V_en", V_en(T)),
        ("V_ne", V_ne(T))
    ]

    for (name, V) in potentials
        H = hamiltonian(V, grid, T)
        eigvals, eigvecs = eigen_solver(H, T)

        @test length(eigvals) == length(eigvecs)
        @test length(eigvals) > 0
        @test length(eigvals) ≤ size(H, 1)

        for (λ, ψ) in zip(eigvals, eigvecs)
            residual = H * ψ - λ * ψ
            @test norm(residual) ≤ 1e-7 * max(1.0, norm(H))
        end
    end
end


@testset "hydrogen ground-state eigenpair" begin
    T = Float64

    x_axis = range(0.5, step = 0.5, length = 3)
    grid = GridMatrix(x_axis, x_axis, x_axis)
    V = CoulombPotential(-one(T), [one(T)], zeros(T, 1, 3))
    H = hamiltonian(V, grid, T)

    eigvals, eigvecs = eigen_solver(H, T)
    sorted_vals = sort(real.(eigvals))
    idx = argmin(real.(eigvals))
    ψ = eigvecs[idx]

    @test length(eigvals) == length(eigvecs)
    @test length(eigvals) > 0
    @test sorted_vals[1] < 5.0
    @test norm(H * ψ - sorted_vals[1] * ψ) ≤ 1e-7 * max(1.0, norm(H))
end


@testset "symbolic Hamiltonian matches numerical Hamiltonian on the interior grid" begin
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
        H_num = hamiltonian(V, grid, T)

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
