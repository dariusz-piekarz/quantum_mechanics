# ============================================================
# Bohr radius
# ============================================================

bohr_radius(::Type{T}) where {T<:AbstractFloat} =
    T(5.2917721054482e-11)

# ============================================================
# radial_function
# ============================================================

"""
    radial_function(n, l, r; factorials=nothing)

Compute the radial part `R_{n,l}(r)` of the hydrogen wave function for array inputs.

This is the broadcasted form used for grid-based evaluations, where `r` can be a
vector, matrix, or higher-dimensional array.
"""
function radial_function(
    n::Integer,
    l::Integer,
    r::AbstractArray{T},
    factorials::Union{Nothing,AbstractVector{<:Integer}}=nothing
)::AbstractArray{<:Number} where {T<:AbstractFloat}

    n < 1 && throw(ArgumentError("n must be >= 1"))
    (l < 0 || l >= n) && throw(ArgumentError("l must be non-negative and < n"))

    fact = resolve_factorials(factorials, n + l)
    α₀ = bohr_radius(T)          # <- dopasowane do typu wejścia, nie na sztywno
    ρ = 1 / (α₀ * n)
    normalizer = sqrt(8 * ρ^3 * fact[n-l] / (2 * n * fact[n+l+1]))

    x = @. 2ρ * r
    L = generalized_laguerre(n - l - 1, 2l + 1, x)
    R_nl = @. normalizer * exp(-x / 2) * x^l * L

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

    α₀ = bohr_radius(Float64)
    fact = resolve_factorials(factorials, n + l)
    ρ = 1 / (α₀ * n)

    R_nl = sqrt(8 * ρ^3 * fact[n-l] / (2 * n * fact[n+l+1])) *
           exp(-ρ * r) * (2 * r * ρ)^l *
           generalized_laguerre(n-l-1, 2l+1, 2 * r * ρ)

    return R_nl
end
