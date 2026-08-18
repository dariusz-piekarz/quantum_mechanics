using Test
using Symbolics
using QuantumMechanics

# ============================================================
# Helpers
# ============================================================

const ATOL = 1e-10
const RTOL = 1e-8

α₀ = bohr_radius(Float64)


# ============================================================
# Radial hydrogen functions
# ============================================================

@testset "Hydrogen radial functions" begin
    r = collect(range(0.0, 20.0 * α₀, length=200))

    for n in 1:4
        for l in 0:(n-1)

            R = radial_function(n, l, r)

            @test size(R) == size(r)
            @test all(isfinite, R)
        end
    end

    # r^l behavior near the origin
    #
    # For l > 0:
    # R_nl(0) = 0
    #
    for n in 2:4
        for l in 1:(n-1)
            R0 = radial_function(n, l, [0.0])
            @test R0[1] ≈ 0.0 atol=ATOL
        end
    end

    # Invalid quantum numbers
    @test_throws ArgumentError radial_function(0, 0, [1.0])
    @test_throws ArgumentError radial_function(2, -1, [1.0])
    @test_throws ArgumentError radial_function(2, 2, [1.0])
end


# ============================================================
# Hydrogen wave function
# ============================================================

@testset "Hydrogen wave function" begin

    r = [0.0, 0.5α₀, α₀, 2α₀]
    θ = [0.0, π/4, π/2, π]
    φ = [0.0, π/2, π, 3π/2]

    ψ = wave_function.(
        2,
        1,
        0,
        r,
        θ,
        φ
    )

    @test size(ψ) == size(r)
    @test eltype(ψ) <: Number
    @test all(isfinite, ψ)

    # Product structure:
    #
    # ψ = R_nl * Y_lm
    #
    facts = factorials_table(4)

    R = radial_function(2, 1, r, facts)
    Y = spherical_harmonic.(1, 0, θ, φ, Ref(facts))

    @test all(isapprox.(ψ, R .* Y; rtol=1e-7, atol=1e-9))

    # Invalid quantum numbers
    @test_throws ArgumentError wave_function(0, 0, 0, r, θ, φ)
    @test_throws ArgumentError wave_function(2, -1, 0, r, θ, φ)
    @test_throws ArgumentError wave_function(2, 2, 0, r, θ, φ)
    @test_throws ArgumentError wave_function(2, 1, 2, r, θ, φ)
end


# ============================================================
# Cartesian wave-function grid
# ============================================================

@testset "Cartesian wave-function grid" begin

    grid = collect(range(-1.0, 1.0, length=5))

    ψ = wave_function_cartesian_grid(
        2,
        1,
        0,
        grid
    )

    @test size(ψ) ==
          (length(grid), length(grid), length(grid))

    @test eltype(ψ) <: Number
    @test all(isfinite, ψ)

    # Check a few points against the spherical-coordinate version.

    test_points = [
        (1, 1, 1),
        (2, 3, 4),
        (3, 2, 1),
    ]

    for (i, j, k) in test_points

        x = grid[i]
        y = grid[j]
        z = grid[k]

        r̃ = sqrt(x^2 + y^2 + z^2)

        if r̃ == 0
            θ = 0.0
            φ = 0.0
        else
            θ = acos(clamp(z / r̃, -1.0, 1.0))
            φ = atan(y, x)
        end

        r = α₀ * r̃

        expected =
            wave_function(
                2,
                1,
                0,
                r,
                θ,
                φ
            )

        @test isapprox(ψ[i, j, k], expected; rtol=1e-6, atol=1e-9)
    end
end


# ============================================================
# Physical normalization of spherical harmonics
# ============================================================

@testset "Spherical harmonic normalization" begin

    # Numerical quadrature using a sufficiently fine tensor grid.
    #
    # ∫ |Y_l^m|² sin(θ) dθ dφ = 1

    θ = collect(range(0.0, π, length=201))
    φ = collect(range(0.0, 2π, length=401))

    dθ = θ[2] - θ[1]
    dφ = φ[2] - φ[1]

    for l in 0:3
        for m in (-l):l

            Y = spherical_harmonic(l, m, θ, φ)

            integral =
                sum(
                    abs2(Y[i, j]) * sin(θ[i])
                    for i in eachindex(θ),
                    j in eachindex(φ)
                ) *
                dθ *
                dφ

            @test integral ≈ 1.0 atol=3e-3
        end
    end
end


@testset "Float32 grid stays Float32 (no silent Float64 promotion)" begin
    grid = collect(range(-5f0, 5f0, length=5))
    @test eltype(wave_function_cartesian_grid(2, 1, 0, grid)) == ComplexF32
    @test eltype(wave_function_cartesian_grid(3, 2, -2, grid)) == ComplexF32
end


@testset "radial_function: eltype stays Float32 for Float32 input" begin
    r = collect(range(0f0, 5f0, length=10))
    @test eltype(radial_function(2, 1, r)) == Float32
    
    r_big = collect(range(0f0, 5f0, length=3))
    @test eltype(radial_function(12, 11, r_big)) isa Type
end
