using Test
using Symbolics

include(joinpath(@__DIR__, "math_special.jl"))
include(joinpath(@__DIR__, "legendre.jl"))
include(joinpath(@__DIR__, "laguerre.jl"))
include(joinpath(@__DIR__, "spherical_harmonics.jl"))
include(joinpath(@__DIR__, "wave_function.jl"))

@testset "factorials and polynomial basics" begin
    @test factorials_table(5) == [1, 1, 2, 6, 24, 120]

    @test legendre(0, 0.5) == 1.0
    @test legendre(1, 0.5) ≈ 0.5
    @test legendre(2, 0.5) ≈ ((3 * 0.5^2 - 1) / 2)
    @test legendre(2, 0.0) ≈ -0.5

    @test laguerre(0, 0.5) == 1.0
    @test laguerre(1, 0.5) ≈ 0.5
    @test laguerre(2, 0.5) ≈ 1.0 - 2 * 0.5 + 0.5^2 / 2
    @test laguerre(2, 2.0) ≈ 1.0 - 4.0 + 2.0^2 / 2.0

    @test generalized_laguerre(0, 0.5, 0.5) == 1.0
    @test generalized_laguerre(1, 0.5, 2.0) ≈ 1.5 - 2.0
end

@testset "spherical harmonics and radial values" begin
    @test spherical_harmonic(0, 0, 0.0, 0.0) ≈ 1 / sqrt(4π)
    @test N(1, 0) ≈ sqrt(3 / (4π))

    radial_r = radial_function(2, 1, [0.1, 0.5])
    @test size(radial_r) == (2,)
    @test all(isfinite, radial_r)

    Y = spherical_harmonic(1, 0, [0.1, 0.2], [0.3, 0.4])
    @test size(Y) == (2, 2)
    @test all(isfinite, Y)
end

@testset "array-based hydrogen wave function" begin
    R = [1.0 2.0; 3.0 4.0]
    θ = fill(0.5, 2, 2)
    φ = fill(0.25, 2, 2)

    ψ = wave_function(3, 2, 0, R, θ, φ)
    @test size(ψ) == size(R)
    @test eltype(ψ) <: Number
    @test all(isfinite, ψ)

    grid = collect(range(-1.0, 1.0, length=3))
    ψ_grid = wave_function_cartesian_grid(3, 1, 0, grid)
    @test size(ψ_grid) == (3, 3, 3)
    @test all(isfinite, ψ_grid)
end
