using Test
using QuantumMechanics
using GLMakie, Meshing, GeometryBasics, LaTeXStrings

@testset "Optional plotting extension" begin
    @test Base.get_extension(QuantumMechanics, :QuantumMechanicsMakieExt) !== nothing
    @test plot_legendre(2; num_points=8) isa GLMakie.Figure
    @test plot_associated_legendre(2, 1; num_points=8) isa GLMakie.Figure
    @test plot_laguerre(2; num_points=8) isa GLMakie.Figure
    @test plot_generalized_laguerre(2, 1; num_points=8) isa GLMakie.Figure
    @test plot_spherical_harmonic(1, 0; num_points=8) isa GLMakie.Figure
end
