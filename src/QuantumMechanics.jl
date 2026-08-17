module QuantumMechanics

using Symbolics
using GLMakie
using LaTeXStrings
using Meshing
using GeometryBasics
using ForwardDiff

# ============================================================
# Source directory
# ============================================================

const SRC_DIR = @__DIR__

# ============================================================
# Mathematical utilities
# ============================================================

include(joinpath(SRC_DIR, "factorials.jl"))

# ============================================================
# Special polynomials
# ============================================================

include(joinpath(SRC_DIR, "legendre.jl"))
include(joinpath(SRC_DIR, "laguerre.jl"))

# ============================================================
# Spherical harmonics
# ============================================================

include(joinpath(SRC_DIR, "spherical_harmonics.jl"))

# ============================================================
# Hydrogen radial function
# ============================================================

include(joinpath(SRC_DIR, "radial_function.jl"))

# ============================================================
# Hydrogen wave functions
# ============================================================

include(joinpath(SRC_DIR, "wave_function.jl"))

# ============================================================
# Isosurfaces
# ============================================================

include(joinpath(SRC_DIR, "isosurface.jl"))

# ============================================================
# Plotting
# ============================================================

include(joinpath(SRC_DIR, "plots.jl"))

# ============================================================
# Converter symbolic functions <-> numeric functions
# ============================================================

include(joinpath(SRC_DIR, "converter.jl"))

# ============================================================
# Benchmarking evaluator for SymbolicFunctions
# ============================================================

include(joinpath(SRC_DIR, "evaluator.jl"))

# ============================================================
# Differential operators
# ============================================================

include(joinpath(SRC_DIR, "differential_operator.jl"))

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

export bohr_radius
export radial_function
export wave_function
export wave_function_cartesian_grid

export show_isosurface
export show_isosurface!

export plot_legendre
export plot_associated_legendre
export plot_spherical_harmonic
export plot_laguerre
export plot_generalized_laguerre
export plot_probability_density
export plot_orbital

export to_function
export SymbolicFunction

export evaluate
export evaluate!

export DifferentialOperator
export D
export gradient
export divergence
export laplacian
export hessian

end