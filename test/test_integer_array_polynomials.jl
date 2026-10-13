using Test
using QuantumMechanics

@testset "Polynomial recurrences with integer arrays" begin
    x = [-1, 0, 1]

    @test legendre(2, x) ≈ [1.0, -0.5, 1.0]
    @test associated_legendre(2, 0, x) ≈ legendre(2, x)
    @test associated_legendre(3, -1, x) ≈
          associated_legendre.(3, -1, x)
    @test generalized_laguerre(2, 0, [0, 1, 2]) ≈ [1.0, -0.5, -1.0]
end

@testset "Radial function preserves Float32" begin
    r = bohr_radius(Float32)
    @test radial_function(2, 1, r) isa Float32
    @test radial_function(2, 1, r) ≈ only(radial_function(2, 1, [r]))
    @test isfinite(radial_function(12, 11, r))
    @test eltype(radial_function(12, 11, [r])) == Float32
    r_high = 20f0 * r
    @test radial_function(20, 19, r_high) > 0f0
    @test only(radial_function(20, 19, [r_high])) > 0f0
    @test radial_function(20, 19, r_high) ≈ only(radial_function(20, 19, [r_high]))
end

@testset "Spherical harmonics promote numeric inputs" begin
    @test spherical_harmonic(1, 0, 0.4f0, 0.2f0) isa ComplexF32
    Y = spherical_harmonic(1, 0, [0, 1], [0, 1])
    @test size(Y) == (2, 2)
    @test eltype(Y) == ComplexF64
    @test all(isfinite, Y)
end
