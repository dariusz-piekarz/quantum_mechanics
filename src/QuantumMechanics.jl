module QuantumMechanics

using Symbolics
using GLMakie
using LaTeXStrings
using Meshing
using GeometryBasics

# ============================================================
# Mathematical utilities
# ============================================================
DIR = dirname(@__DIR__)

include(joinpath(DIR, "..", "math_special.jl"))

# ============================================================
# Special polynomials
# ============================================================

include(joinpath(DIR, "..", "legendre.jl"))
include(joinpath(DIR, "..", "laguerre.jl"))

# ============================================================
# Spherical harmonics
# ============================================================

include(joinpath(DIR, "..", "spherical_harmonics.jl"))

# ============================================================
# Hydrogen wave functions
# ============================================================

include(joinpath(DIR, "..", "wave_function.jl"))

# ============================================================
# Isosurfaces
# ============================================================

include(joinpath(DIR, "..", "isosurface.jl"))

# ============================================================
# Plotting
# ============================================================

include(joinpath(DIR, "..", "plots.jl"))


# ============================================================
# Exports
# ============================================================

export factorials_table
export resolve_factorials

export legendre
export associated_legendre

export laguerre
export generalized_laguerre

export N
export P
export spherical_harmonic

export radial_function
export wave_function
export wave_function_cartesian_grid

export show_isosurface
export show_isosurface!

export plot_legendre
export plot_associated_legendre
export plot_spherical_harmonic
export plot_laguerre
export plot_associated_laguerre

export plot_wave_function
export plot_probability_density
export plot_orbital

end