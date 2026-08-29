module QuantumMechanics

using Symbolics
using GLMakie
using LaTeXStrings
using Meshing
using GeometryBasics
using ForwardDiff
using SparseArrays
using KrylovKit

# ============================================================
# Source directory
# ============================================================

const SRC_DIR = @__DIR__
const SPEC_FUNCT_DIR = joinpath(SRC_DIR, "special_functions")
const HYDROGEN_ORB_DIR = joinpath(SRC_DIR, "hydrogen_orbitals")
const VISUALISATION_DIR = joinpath(SRC_DIR, "visualisation")
const SYMB_FUNCT_DIR = joinpath(SRC_DIR, "symbolic_functions")
const SYMB_OP_DIR = joinpath(SRC_DIR, "symbolic_operators")
const NUM_SOLV_DIR = joinpath(SRC_DIR, "numerical_solver")

# ============================================================
# Mathematical utilities
# ============================================================

include(joinpath(SPEC_FUNCT_DIR, "factorials.jl"))

# ============================================================
# Special polynomials
# ============================================================

include(joinpath(SPEC_FUNCT_DIR, "legendre.jl"))
include(joinpath(SPEC_FUNCT_DIR, "laguerre.jl"))

# ============================================================
# Spherical harmonics
# ============================================================

include(joinpath(SPEC_FUNCT_DIR, "spherical_harmonics.jl"))

# ============================================================
# Hydrogen radial function
# ============================================================

include(joinpath(HYDROGEN_ORB_DIR, "radial_function.jl"))

# ============================================================
# Hydrogen wave functions
# ============================================================

include(joinpath(HYDROGEN_ORB_DIR, "wave_function.jl"))

# ============================================================
# Isosurfaces
# ============================================================

include(joinpath(VISUALISATION_DIR, "isosurface.jl"))

# ============================================================
# Plotting
# ============================================================

include(joinpath(VISUALISATION_DIR, "plots.jl"))

# ============================================================
# Converter symbolic functions <-> numeric functions
# ============================================================

include(joinpath(SYMB_FUNCT_DIR, "converter.jl"))

# ============================================================
# Benchmarking evaluator for SymbolicFunctions
# ============================================================

include(joinpath(SYMB_FUNCT_DIR, "evaluator.jl"))

# ============================================================
# Differential operators
# ============================================================

include(joinpath(SYMB_OP_DIR, "differential_operator.jl"))

# ============================================================
# Potential
# ============================================================

include(joinpath(SYMB_OP_DIR, "potential.jl"))

# ============================================================
# Kinetic Operator
# ============================================================

include(joinpath(SYMB_OP_DIR, "kinetic_operator.jl"))

# ============================================================
# Hamiltonian
# ============================================================

include(joinpath(SYMB_OP_DIR, "hamiltonian.jl"))

# ============================================================
# discretized laplacian
# ============================================================

include(joinpath(NUM_SOLV_DIR, "discretized_laplacian.jl"))

# ============================================================
# eigensolver
# ============================================================

include(joinpath(NUM_SOLV_DIR, "eigensolver.jl"))

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

export AbstractPotential
export CoulombPotential
export HarmonicPotential
export V
export V_e
export V_ee
export V_en
export V_ne

export AbstractOperator
export KineticOperator
export T̂

export Hamiltonian

export raw_laplacian
export eigen_solver

end