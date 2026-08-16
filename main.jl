DIR = joinpath(@__DIR__, "src")

include(joinpath(DIR, "special_functions", "math_special.jl"))
include(joinpath(DIR, "special_functions", "legendre.jl"))
include(joinpath(DIR, "special_functions", "laguerre.jl"))
include(joinpath(DIR, "special_functions", "spherical_harmonics.jl"))
include(joinpath(DIR, "hydrogen", "radial_function.jl"))
include(joinpath(DIR, "hydrogen", "wave_function.jl"))
include(joinpath(DIR, "visualisation", "plots.jl"))

#---------------------PARAMS---------------------
n, l, m = 3, 2, -2
probability_lvl = Float32(0.1)
@show "--- New run: n: $n, l: $l, m: $m ---"

#---------------------SYMBOLICS---------------------
@variables x, α, r, θ, φ
#@time legendre(2, x)
#@time associated_legendre(4, 2, x)
#@time laguerre.(4, x)
#@time generalized_laguerre(2, α, x)
#@time generalized_laguerre(2, 0.4, x)
#@time spherical_harmonic(3, 2, θ, φ)
#@time wave_function(n, m, l, r, θ, φ)

#---------------------NUMERIC---------------------
#@time legendre(2, 0.2)
#@time generalized_laguerre(2, 0.4, 2)
#@time spherical_harmonic(3, 2, 0.2, 3.1)
#@time wave_function(n, m, l, 0.1, 0.34, 0.66)

#---------------------VECTOR---------------------
#@time associated_legendre(4, 2, [0.5, 0.55])
#@time laguerre.(4, [7,8])
#@time spherical_harmonic.(3, 2, [0.2, 0.4, 0.4], [3.1, 1, 2])
#@time wave_function(n, m, l, [0.2, 0.4], [1, 2], [2.4, 2.9, 4.5, 3.7])

#---------------------PLOTS---------------------
#plot_legendre(5)
#plot_associated_legendre(4, 2)
#plot_generalized_laguerre(9, 1.5)
#plot_spherical_harmonic(8, 7, draw=:real)
#plot_probability_density(n, l, m, level=probability_lvl)
plot_orbital(n, l, m, level=probability_lvl)