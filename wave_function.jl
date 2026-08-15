include(joinpath(@__DIR__, "spherical_harmonics.jl"))
include(joinpath(@__DIR__, "laguerre.jl"))

α₀ = 5.2917721054482*10^(-11)

"""
    radial_function(n, l, r; factorials=nothing)

Compute the radial part `R_{n,l}(r)` of the hydrogen wave function for array inputs.

This is the broadcasted form used for grid-based evaluations, where `r` can be a
vector, matrix, or higher-dimensional array.
"""
function radial_function(
    n::Integer,
    l::Integer,
    r::AbstractArray{<:Number},
    factorials::Union{Nothing,AbstractVector{<:Integer}}=nothing
)::AbstractArray{<:Number}
    n < 1 && throw(ArgumentError("n must be >= 1"))
    (l < 0 || l >= n) && throw(ArgumentError("l must be non-negative and < n"))

    fact = resolve_factorials(factorials, n + l)
    ρ = 1 / (α₀ * n)
    normalizer = sqrt(8 * ρ^3 * fact[n-l] / (2 * n * fact[n+l+1]))

    R_nl = normalizer .* exp.(-ρ .* r) .* (2 * ρ .* r) .^ l .* generalized_laguerre(n-l-1, 2l+1, 2 * ρ .* r)
    return R_nl
end


"""
    radial_function(n, l, r; factorials=nothing)

Compute the radial part `R_{n,l}(r)` for a single radial coordinate.
"""
function radial_function(
    n::Integer,
    l::Integer,
    r::Union{Num,Number},
    factorials::Union{Nothing,AbstractVector{<:Integer}}=nothing
)::Union{Num,Number}

    n < 1 && throw(ArgumentError("n must be >= 1"))
    (l < 0 || l >= n) && throw(ArgumentError("l must be non-negative and < n"))

    fact = resolve_factorials(factorials, n + l)
    ρ = 1 / (α₀ * n)

    R_nl = sqrt(8 * ρ^3 * fact[n-l] / (2 * n * fact[n+l+1])) *
           exp(-ρ * r) * (2 * r * ρ)^l *
           generalized_laguerre(n-l-1, 2l+1, 2 * r * ρ)
    return R_nl
end




"""
    wave_function(n, l, m, r, θ, φ; factorials=nothing)

Evaluate the hydrogen wave function on a full array-based grid.

The result is the elementwise product of the radial term and the spherical
harmonic, with broadcasting over all grid points.
"""
function wave_function(
    n::Integer,
    l::Integer,
    m::Integer,
    r::AbstractArray{<:Number},
    θ::AbstractArray{<:Number},
    φ::AbstractArray{<:Number},
    factorials::Union{Nothing,AbstractVector{<:Integer}}=nothing
)::AbstractArray{<:Number}

    facts = resolve_factorials(factorials, max(2l, n+l))
    R_nl = radial_function(n, l, r, facts)
    Y_lm = spherical_harmonic(l, m, θ, φ, facts)
    return R_nl .* Y_lm
end


"""
    wave_function(n::Integer, l::Integer, m::Integer, r::Union{Num,Number}, θ::Union{Num,Number}, φ::Union{Num,Number})

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
    φ::Union{Num,Number},
    factorials::Union{Nothing,AbstractVector{<:Integer}}=nothing
)::Union{Num,Number}
    facts = resolve_factorials(factorials, max(2l, n+l))
    R_nl = radial_function(n, l, r, facts)

    Y_lm = spherical_harmonic(l, m, θ, φ, facts)
    return R_nl * Y_lm
end


"""
    wave_function_cartesian_grid(n, l, m, ξ; factorials=nothing)

Build a Cartesian grid from `ξ`, compute spherical coordinates, and evaluate the
hydrogen wave function on the resulting 3D mesh.
"""
function wave_function_cartesian_grid(
    n::Integer,
    l::Integer,
    m::Integer,
    ξ::AbstractVector;
    factorials::Union{Nothing,AbstractVector{<:Integer}}=nothing
)

    ξ = collect(Float64, ξ)
    N = length(ξ)

    x = reshape(ξ, :, 1, 1)
    y = reshape(ξ, 1, :, 1)
    z = reshape(ξ, 1, 1, :)

    X = repeat(x, 1, N, N)
    Y = repeat(y, N, 1, N)
    Z = repeat(z, N, N, 1)

    R2 = X .^ 2 .+ Y .^ 2 .+ Z .^ 2
    R = α₀ .* sqrt.(R2)
    denom = sqrt.(R2) .+ eps(Float64)
    θ = @. ifelse(R == 0, 0.0, acos(clamp(Z / denom, -1.0, 1.0)))
    φ = @. ifelse(R == 0, 0.0, atan(Y, X))

    fact = resolve_factorials(factorials, max(2l, n+l))
    return wave_function(n, l, m, R, θ, φ, fact)
end
