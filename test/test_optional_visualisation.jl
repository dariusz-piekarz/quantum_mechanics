using Test
using QuantumMechanics

@testset "Plotting is optional" begin
    @test Base.get_extension(QuantumMechanics, :QuantumMechanicsMakieExt) === nothing
    @test_throws ArgumentError plot_legendre(2)
end
