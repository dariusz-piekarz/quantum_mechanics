# QuantumMechanics.jl

Eksperymentalny pakiet Julia do obliczeń związanych z mechaniką kwantową. Zawiera wielomiany specjalne, harmoniki sferyczne i orbitale atomu wodoru, operatory symboliczne, rzadkie dyskretyzacje hamiltonianu oraz wizualizacje oparte na GLMakie.

## Uruchomienie

Wymagana jest Julia 1.12. W katalogu repozytorium:

```sh
julia --project=. -e "using Pkg; Pkg.instantiate()"
julia --project=. -e "using QuantumMechanics; println(legendre(2, 0.5))"
```

Przykłady znajdują się w katalogu `example/`. Na początek można uruchomić:

```sh
julia --project=. example/example_symbolic_calculations.jl
julia --project=. example/example_numerical_calculations.jl
```

Przykłady w `example/` z `visualisations` w nazwie otwierają okna graficzne.

## Testy

```sh
julia --project=. -e "using Pkg; Pkg.test()"
```

## Jednostki

`bohr_radius(Float64)` zwraca promień Bohra w metrach; `radial_function` i `wave_function` przyjmują odległość w metrach. `wave_function_cartesian_grid` przyjmuje wektor współrzędnych w wielokrotnościach promienia Bohra. Dyskretny hamiltonian ma postać `-½Δ + V` i nie wykonuje przeliczania jednostek. Przy użyciu potencjału Coulomba `V = -1/r` należy podawać współrzędne siatki w jednostkach atomowych, aby energie były w hartree.

## Granice i kroki siatki numerycznej

`GridMatrix(min, max, h)` tworzy siatkę sześcienną. `min` jest pierwszą próbką, a `max` górną granicą **pozycji próbek**: ostatnia próbka nigdy nie przekracza `max`. Jeżeli `(max-min)/h` nie jest całkowite, `max` nie należy do siatki. Przykładowo `GridMatrix(0.0, 1.0, 0.4)` ma na każdej osi punkty `0.0, 0.4, 0.8`. Każda oś musi zawierać co najmniej dwa punkty.

Różne kroki i liczby punktów osi można podać jako `GridMatrix(x_range, y_range, z_range)` albo `GridMatrix(xmin, xmax, ymin, ymax, zmin, zmax, hx, hy, hz)`. W `grid.grid` kolumny oznaczają `x, y, z`, a `x` zmienia się najszybciej. Współrzędne, granice i kroki muszą mieć zgodne jednostki. Dla `hamiltonian` są to bohr, a otrzymane energie są w hartree, o ile również potencjał jest wyrażony w jednostkach atomowych.

Operator różnic skończonych zakłada zerową funkcję falową **o jeden krok poza skrajnymi próbkami** na każdej osi. Dla stałego fizycznego pudełka z warunkiem Dirichleta na `0` i `L` należy więc próbkować wnętrze: `h = L/(n+1)` i `axis = range(h; step=h, length=n)`. Przy zagęszczaniu siatki zmieniamy `n` i `h`, ale utrzymujemy te same fizyczne granice pudełka.

## Zakres API

- Funkcje specjalne: `legendre`, `associated_legendre`, `laguerre`, `generalized_laguerre`, `spherical_harmonic`.
- Wodór: `radial_function`, `wave_function`, `wave_function_cartesian_grid`.
- Symbole i operatory: `SymbolicFunction`, `gradient`, `laplacian`, `KineticOperator`, `Hamiltonian`.
- Obliczenia na siatce: `GridMatrix`, `hamiltonian`, `eigen_solver`.

Szczegóły wywołań pokazują skrypty w `example/` oraz testy w `test/`.
