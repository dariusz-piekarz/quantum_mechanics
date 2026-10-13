module QuantumMechanicsMakieExt

using QuantumMechanics: factorials_table, legendre, associated_legendre,
    spherical_harmonic, laguerre, generalized_laguerre, wave_function_cartesian_grid
import QuantumMechanics: show_isosurface, show_isosurface!, plot_legendre,
    plot_associated_legendre, plot_spherical_harmonic, plot_laguerre,
    plot_generalized_laguerre, plot_probability_density, plot_orbital

using GLMakie
using LaTeXStrings
using Meshing
using GeometryBasics

include("visualisation/isosurface.jl")
include("visualisation/plots.jl")

end
