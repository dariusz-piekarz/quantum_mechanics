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

include("special_functions/factorials.jl")

# ============================================================
# Special polynomials
# ============================================================

include("special_functions/legendre.jl")
include("special_functions/laguerre.jl")

# ============================================================
# Spherical harmonics
# ============================================================

include("special_functions/spherical_harmonics.jl")

# ============================================================
# Hydrogen radial function
# ============================================================

include("hydrogen_orbitals/radial_function.jl")

# ============================================================
# Hydrogen wave functions
# ============================================================

include("hydrogen_orbitals/wave_function.jl")

# ============================================================
# Isosurfaces
# ============================================================

include("visualisation/isosurface.jl")

# ============================================================
# Plotting
# ============================================================

include("visualisation/plots.jl")

# ============================================================
# Converter symbolic functions <-> numeric functions
# ============================================================

include("symbolic_functions/converter.jl")

# ============================================================
# Benchmarking evaluator for SymbolicFunctions
# ============================================================

include("symbolic_functions/evaluator.jl")

# ============================================================
# Differential operators
# ============================================================

include("symbolic_operators/differential_operator.jl")

# ============================================================
# Potential
# ============================================================

include("symbolic_operators/potential.jl")

# ============================================================
# Kinetic Operator
# ============================================================

include("symbolic_operators/kinetic_operator.jl")

# ============================================================
# Hamiltonian
# ============================================================

include("symbolic_operators/hamiltonian.jl")

# ============================================================
# discretized laplacian
# ============================================================

include("numerical_solver/discretized_laplacian.jl")

# ============================================================
# discretized hamiltonian
# ============================================================

include("numerical_solver/discretized_hamiltonian.jl")

# ============================================================
# eigensolver
# ============================================================

include("numerical_solver/eigensolver.jl")

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

export Grid
export GridMatrix
export grid_step
export grid_size
export raw_laplacian_1d
export laplacian_1d
export kinetic_operator
export potential_operator
export hamiltonian
export eigen_solver

end