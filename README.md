# QuantumMechanics.jl

An experimental Julia package for quantum mechanics calculations. It includes special polynomials, spherical harmonics and hydrogen atom orbitals, symbolic operators, sparse discretizations of the Hamiltonian, and optional visualizations built with GLMakie.

## Getting started

Julia 1.12 is required. Run these commands from the repository root:

```sh
julia --project=. -e "using Pkg; Pkg.instantiate()"
julia --project=. -e "using QuantumMechanics; println(legendre(2, 0.5))"
```

Examples are in `example/`. These are good starting points:

```sh
julia --project=. example/example_symbolic_calculations.jl
julia --project=. example/example_numerical_calculations.jl
```

The numerical examples and core tests do not require GLMakie or a display. Plotting lives in the optional `QuantumMechanicsMakieExt` extension. To run the graphical examples, instantiate the separate environment and select it:

```sh
julia --project=visualisations -e "using Pkg; Pkg.instantiate()"
julia --project=visualisations example/example_visualisations.jl
julia --project=visualisations example/example_hydrogen_atom_orbital_visualisations.jl
```

The graphical scripts load all four optional packages (`GLMakie`, `Meshing`, `GeometryBasics`, `LaTeXStrings`), which activates the extension. If using another environment, add and load those packages there before calling `plot_*` or `show_isosurface*`.

## Tests

```sh
julia --project=. -e "using Pkg; Pkg.test()"
```

An optional plotting smoke test is available after instantiating `visualisations/`:

```sh
julia --project=visualisations visualisations/smoke_test.jl
```

## Units

`bohr_radius(Float64)` returns the Bohr radius in meters; `radial_function` and `wave_function` take distances in meters. `wave_function_cartesian_grid` takes coordinates expressed as multiples of the Bohr radius. The discretized Hamiltonian has the form `-½Δ + V` and performs no unit conversion. For a Coulomb potential `V = -1/r`, supply grid coordinates in atomic units (bohr) to obtain energies in hartree.

## Numerical grid bounds and spacing

`GridMatrix(min, max, h)` creates a cubic grid. `min` is the first sample, while `max` is an inclusive upper bound on **sample positions**: the last sample never exceeds `max`. If `(max-min)/h` is not an integer, `max` is not a grid point. For example, `GridMatrix(0.0, 1.0, 0.4)` gives `0.0, 0.4, 0.8` on each axis. Each axis must contain at least two points.

For different step sizes or point counts along the axes, use `GridMatrix(x_range, y_range, z_range)` or `GridMatrix(xmin, xmax, ymin, ymax, zmin, zmax, hx, hy, hz)`. The columns of `grid.grid` are `x, y, z`, with `x` varying fastest. Coordinates, bounds, and steps must use consistent units. For `hamiltonian`, coordinates should be in bohr; energies are then in hartree when the potential is also expressed in atomic units.

The finite difference operator assumes a zero wave function **one step beyond the outermost samples** on each axis. To keep a physical box with Dirichlet boundaries at `0` and `L` fixed while refining the grid, sample its interior with `h = L/(n+1)` and `axis = range(h; step=h, length=n)`. Change `n` and `h` while keeping the physical box boundaries fixed.

## API overview

- Special functions: `legendre`, `associated_legendre`, `laguerre`, `generalized_laguerre`, `spherical_harmonic`.
- Hydrogen atom: `radial_function`, `wave_function`, `wave_function_cartesian_grid`.
- Symbolic functions and operators: `SymbolicFunction`, `gradient`, `laplacian`, `KineticOperator`, `Hamiltonian`.
- Grid calculations: `GridMatrix`, `hamiltonian`, `eigen_solver`.

The scripts in `example/` and the tests in `test/` demonstrate their usage.

## Development roadmap

The detailed solver and visualization roadmap is in [to_do.md](to_do.md) (written in Polish).
