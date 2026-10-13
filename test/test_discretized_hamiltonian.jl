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

    shifted = GridMatrix(
        x,
        range(0.3, step = 0.1, length = 3),
        range(1.0, step = 0.1, length = 3),
    )
    @test size(shifted.grid) == (27, 3)
    @test size(kinetic_operator(shifted, Float64)) == (27, 27)

    @test_throws DimensionMismatch GridMatrix(range(0.0, step = 0.1, length = 1), y, z)
    different_steps = GridMatrix(x, range(0.0, step = 0.2, length = 3), z)
    @test size(kinetic_operator(different_steps, T)) == (27, 27)
    @test_throws DimensionMismatch GridMatrix{Float64}(ones(8, 2))

    # The point count alone cannot reveal the axis dimensions.
    rectangular = GridMatrix(
        range(0.0, step = 0.5, length = 2),
        range(0.0, step = 0.5, length = 4),
        range(0.0, step = 0.5, length = 8),
    )
    @test size(kinetic_operator(rectangular, Float64)) == (64, 64)

    bounded = GridMatrix(0.0, 1.0, 0.4)
    @test size(bounded.grid) == (27, 3)
    @test bounded.grid[end, :] ≈ [0.8, 0.8, 0.8]
    @test all(bounded.grid .<= 1.0)
    @test GridMatrix(0.0, 0.3, 0.1).grid[end, 1] ≈ 0.3
    @test GridMatrix(0, 2, 1) isa GridMatrix{Float64}

    anisotropic = GridMatrix(0.0, 0.5, 0.0, 0.5, 0.0, 1.0, 0.5, 0.25, 1.0)
    @test size(anisotropic.grid) == (12, 3)
    K = kinetic_operator(anisotropic, T)
    @test K[1, 1] ≈ 21.0
    @test K[1, 2] ≈ -2.0  # x neighbor
    @test K[1, 3] ≈ -8.0  # y neighbor
    @test K[1, 7] ≈ -0.5  # z neighbor
    @test K == transpose(K)

    irregular = GridMatrix{Float64}(copy(anisotropic.grid))
    irregular.grid[5, 2] += 0.01
    @test_throws DimensionMismatch kinetic_operator(irregular, T)

    @test_throws DimensionMismatch GridMatrix(0.0, 0.1, 0.2)
    @test_throws DimensionMismatch GridMatrix(0.0, 1.0, 0.0)
    @test_throws ArgumentError GridMatrix(0.0, 1.0, Inf)
end


@testset "free-particle energy converges as grid spacing decreases" begin
    # The sampled points are strictly inside the fixed Dirichlet box [0, L].
    # The zero-valued ghost points lie exactly at 0 and L on every axis.
    L = 4.0
    exact_energy = 3π^2 / (2L^2)
    errors = Float64[]
    for n in 2:5
        h = L / (n + 1)
        grid = GridMatrix(range(h; step=h, length=n))
        energy = first(eigen_solver(kinetic_operator(grid, Float64), Float64; howmany=1)[1])
        push!(errors, abs(energy - exact_energy))
    end
    @test all(diff(errors) .< 0)
    @test last(errors) < first(errors) / 3

    # Rectangular boxes retain fixed boundaries while each axis has its own h.
    lengths = (4.0, 5.0, 6.0)
    exact_rectangular = π^2 / 2 * sum(inv(Li^2) for Li in lengths)
    rectangular_errors = Float64[]
    for dims in ((2, 3, 2), (3, 4, 3))
        axes = ntuple(3) do i
            h = lengths[i] / (dims[i] + 1)
            range(h; step=h, length=dims[i])
        end
        grid = GridMatrix(axes...)
        energy = first(eigen_solver(kinetic_operator(grid, Float64), Float64; howmany=1)[1])
        push!(rectangular_errors, abs(energy - exact_rectangular))
    end
    @test rectangular_errors[2] < rectangular_errors[1]
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


@testset "eigen_solver resolves low and degenerate eigenvalues" begin
    A = spdiagm(0 => [3.0, 2.0, 1.0])
    vals, vecs = eigen_solver(A, Float64; howmany=2)
    @test sort(vals) ≈ [1.0, 2.0]
    @test all(norm(A * vec - val * vec) < 1e-8 for (val, vec) in zip(vals, vecs))
    high_vals, _ = eigen_solver(A, Float64; howmany=2, which=:LR)
    @test high_vals ≈ [3.0, 2.0]

    repeated = spdiagm(0 => [1.0, 1.0, 3.0])
    repeated_vals, repeated_vecs = eigen_solver(repeated, Float64; howmany=2)
    @test sort(repeated_vals) ≈ [1.0, 1.0]
    @test length(repeated_vecs) == 2

    @test_throws ArgumentError eigen_solver(A, Float64; howmany=0)
    @test_throws DimensionMismatch eigen_solver(ones(2, 3), Float64)
    @test_throws ArgumentError eigen_solver([1.0 2.0; 3.0 4.0], Float64)
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
    grid_shape = (length(x_axis), length(x_axis), length(x_axis))

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

        H_sym_flat = Vector{T}(undef, size(grid.grid, 1))
        evaluate!(H_sym, grid.grid, H_sym_flat)

        # The discrete finite-difference Laplacian is only expected to match the
        # symbolic operator on the interior of the grid. The boundary stencil is
        # a truncated approximation and therefore differs by construction.
        H_num_values = reshape(H_num * ψ_values, grid_shape)
        H_sym_values = reshape(H_sym_flat, grid_shape)

        H_num_interior = H_num_values[2:end-1, 2:end-1, 2:end-1]
        H_sym_interior = H_sym_values[2:end-1, 2:end-1, 2:end-1]

        @test size(H_num_interior) == (1, 1, 1)
        @test H_num_interior ≈ H_sym_interior atol = 1e-10 rtol = 1e-10
    end
end
