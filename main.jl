include(joinpath(@__DIR__, "plots.jl"))


@variables x, α
#associated_legendre(2,0,[1.5, 2.0, 2.5])
#generalized_lagguerre(2, α, x)
plot_generalized_laguerre(2, 0.5)
#plot_spherical_harmonic(6, 0, draw=:probability)
#plot_legendre(5)
#plot_associated_legendre(5, 3)

