# ============================================================
# Bohr radius
# ============================================================

bohr_radius(::Type{T}) where {T<:AbstractFloat} =
    T(5.2917721054482e-11)

# Keep the factorial ratio in a wider type until after the square root. For
# large quantum numbers the ratio alone can underflow in Float32 even though
# the complete normalization is representable.
function _radial_normalizer(n::Integer, l::Integer, ρ::T, fact::AbstractVector{<:Integer}) where {T<:AbstractFloat}
    W = n + l > 20 ? BigFloat : promote_type(T, Float64)
    ρ_wide = W(ρ)
    numerator = W(8) * ρ_wide^3 * W(fact[n-l])
    denominator = W(2n) * W(fact[n+l+1])
    return T(sqrt(numerator / denominator))
end

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
    α₀ = bohr_radius(T)       
    ρ = 1 / (α₀ * n)
    normalizer = _radial_normalizer(n, l, ρ, fact)

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

    T = r isa AbstractFloat ? typeof(r) : Float64
    α₀ = bohr_radius(T)
    fact = resolve_factorials(factorials, n + l)
    ρ = 1 / (α₀ * n)

    normalizer = _radial_normalizer(n, l, ρ, fact)
    R_nl = normalizer *
           exp(-ρ * r) * (2 * r * ρ)^l *
           generalized_laguerre(n-l-1, 2l+1, 2 * r * ρ)

    return R_nl
end
