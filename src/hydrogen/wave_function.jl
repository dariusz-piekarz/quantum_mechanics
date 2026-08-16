# ============================================================
# the wave function
# ============================================================

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

    t1 = time()
    R_nl = radial_function(n, l, r, facts)
    t2 = time()
    @info "radial_function: $(t2-t1)"

    t1 = time()
    Y_lm = spherical_harmonic(l, m, θ, φ, facts)
    t2 = time()
    @info "spherical harmonics: $(t2-t1)"

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

# ============================================================
# the wave fucntion cartesian grid
# ============================================================

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
    α₀ = bohr_radius(Float32)
    ξ = collect(Float32, ξ)

    x = reshape(ξ, :, 1, 1)
    y = reshape(ξ, 1, :, 1)
    z = reshape(ξ, 1, 1, :)

    R2 = @. x^2 + y^2 + z^2
    R = @. α₀ * sqrt(R2)
    denom = @. sqrt(R2) + eps(Float32)
    θ = @. ifelse(R2 == 0f0, 0f0, acos(clamp(z / denom, -1f0, 1f0)))
    φ = @. ifelse(R2 == 0f0, 0f0, atan(y, x))
    fact = resolve_factorials(factorials, max(2l, n+l))

    return wave_function(n, l, m, R, θ, φ, fact)
end
