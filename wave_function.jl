include(joinpath(@__DIR__, "spherical_harmonics.jl"))
include(joinpath(@__DIR__, "laguerre.jl"))

α₀ = 5.2917721054482*10^(-11)

"""
    wave_function(n::Integer, l::Integer, m::Integer, r, θ, φ)

Calculate the hydrogen atom wave function (eigenfunction of the Schrödinger equation).

The wave function is a product of radial and angular parts:
ψ(r, θ, φ) = R_{n,l}(r) · Y_l^m(θ, φ)

# Arguments
- `n::Integer`: Principal quantum number (n ≥ 1)
- `l::Integer`: Orbital angular momentum quantum number (0 ≤ l < n)
- `m::Integer`: Magnetic quantum number (-l ≤ m ≤ l)
- `r::Union{Num,Number}`: Radial distance from nucleus
- `θ::Union{Num,Number}`: Polar angle (colatitude)
- `φ::Union{Num,Number}`: Azimuthal angle

# Returns
- Wave function value ψ(r, θ, φ)

# Throws
- `ArgumentError`: If quantum numbers violate selection rules

# References
Bohr radius: α₀ ≈ 5.29 × 10⁻¹¹ m (0.53 Å)
"""
function wave_function(
    n::Integer,
    l::Integer,
    m::Integer,
    r::Union{Num,Number},
    θ::Union{Num,Number},
    φ::Union{Num,Number}
)::Union{Num,Number}
    n < 1 && throw(ArgumentError("n must be >= 1"))
    (l < 0 || l >= n) && throw(ArgumentError("l must be non-negative and < n"))
    abs(m) > l && throw(ArgumentError("|m| must be <= l"))

    ρ = 1 / (α₀ * n)

    R_nl = sqrt(8 * ρ^3 * factorial(n-l-1) / (2 * n * factorial(n+l))) *
           exp(-ρ*r) * (2*r*ρ)^l *
           generalized_lagguerre(n-l-1, 2l+1, 2*r*ρ)

    Y_lm = spherical_harmonic(l, m, θ, φ)

    return R_nl * Y_lm
end


"""
    wave_function(n::Integer, l::Integer, m::Integer, r::AbstractVector, θ::AbstractVector, φ::AbstractVector)

Calculate the hydrogen atom wave function on vectorized spatial grids.

Computes ψ(r, θ, φ) for multiple radial and angular coordinates, returning a 3D array
with dimensions (length(r), length(θ), length(φ)).

# Arguments
- `n::Integer`: Principal quantum number (n ≥ 1)
- `l::Integer`: Orbital angular momentum quantum number (0 ≤ l < n)
- `m::Integer`: Magnetic quantum number (-l ≤ m ≤ l)
- `r::AbstractVector{<:Number}`: Vector of radial distances from nucleus
- `θ::AbstractVector{<:Number}`: Vector of polar angles (colatitude)
- `φ::AbstractVector{<:Number}`: Vector of azimuthal angles

# Returns
- 3D array of shape (length(r), length(θ), length(φ)) containing wave function values

# Example
```julia
r = range(0.1, 10.0, length=50)
θ = range(0, π, length=40)
φ = range(0, 2π, length=60)
ψ = wave_function(2, 1, 0, r, θ, φ)
```
"""
function wave_function(
    n::Integer,
    l::Integer,
    m::Integer,
    r::AbstractVector{<:Number},
    θ::AbstractVector{<:Number},
    φ::AbstractVector{<:Number},
)

    # Radial scaling factor
    ρ = 1 / (α₀ * n)

    # Radial part: R_{n,l}(r) - applies to all r values
    R = @. sqrt(8 * ρ^3 * factorial(n-l-1) / (2 * n * factorial(n+l))) *
           exp(-ρ*r) * (2*r*ρ)^l *
           generalized_lagguerre(n-l-1, 2l+1, 2*r*ρ)

    # Angular part: Y_l^m(θ, φ) - applies to all θ, φ values
    Y = spherical_harmonic(l, m, θ, φ)

    # Full wave function: ψ(r,θ,φ) = R_{n,l}(r) · Y_l^m(θ, φ)
    # Broadcasting to create 3D array
    return reshape(R, :, 1, 1) .*
           reshape(Y, 1, size(Y, 1), size(Y, 2))
end