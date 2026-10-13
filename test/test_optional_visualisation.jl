using Test
using QuantumMechanics

@testset "Plotting is optional" begin

    ext = Base.get_extension(
        QuantumMechanics,
        :QuantumMechanicsMakieExt
    )

    if isnothing(ext)
        # Extension was not loaded, so the plotting function should throw an error.
        @test_throws ArgumentError plot_legendre(2)
    else
        # Extension was loaded, so the plotting function should work.
        @test ext isa Module
        @test hasmethod(plot_legendre, Tuple{Int})
    end

end
