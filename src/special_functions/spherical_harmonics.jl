using Logging

# ============================================================
# Spherical harmonic Normalization const
# ============================================================

"""
    N(l, m; factorials=nothing)

Return the normalization constant for the spherical harmonic `Y_l^m`.

For the standard complex convention,

N_{l,m} = sqrt((2l + 1) / (4π) * (l - |m|)! / (l + |m|)!)
for `m >= 0`, while for negative `m` the same signed convention is used with the
usual relation `Y_l^{-m} = (-1)^m conj(Y_l^m)`, which requires

N_{l,-m} = sqrt((2l + 1) / (4π) * (l + |m|)! / (l - |m|)!).

If `factorials` is provided, it must be a factorial table with convention
`factorials[i+1] = i!`.
"""
function N(
    l::Integer,
    m::Integer,
    factorials::Union{Nothing,AbstractVector{<:Integer}}=nothing
)::Number

    l < 0 && throw(ArgumentError("l must be non-negative"))
    abs(m) > l && throw(ArgumentError("|m| must be <= l"))

    m_abs = abs(m)
    fact = resolve_factorials(factorials, l + m_abs)

    idx_left = l - m_abs + 1
    idx_right = l + m_abs + 1
    if idx_left > length(fact) || idx_right > length(fact)
        throw(ArgumentError("factorials vector is too short for the requested (l, m)"))
    end

    l_m_mabs = fact[idx_left]
    l_p_m_abs = fact[idx_right]

    if m >= 0
        return sqrt((2l + 1) / (4π) * l_m_mabs / l_p_m_abs)
    end

    return sqrt((2l + 1) / (4π) * l_p_m_abs / l_m_mabs)
end

# ============================================================
# the azimuthal phase P
# ============================================================

"""
    P(m, φ)

Return the azimuthal phase factor `exp(i*m*φ)` for an array of angles.
"""
function P(m::Integer, φ::AbstractArray{<:Number})::AbstractArray{<:Number}
    return exp.(im .* m .* φ)
end


"""
    P(m, φ)

Return the azimuthal phase factor `exp(i*m*φ)` for the value in `φ`.
"""
function P(m::Integer, φ::Number)::Number
    return exp(im * m * φ)
end


"""
    P(m, φ)

Return the azimuthal phase factor `exp(i*m*φ)` for the value in `φ`.
"""
function P(m::Integer, φ::Num)::Num
    return exp(im * m * φ)
end


# ============================================================
# Spherical harmonics Y_l^m(θ, φ)
# ============================================================

"""
    spherical_harmonic(l, m, θ, φ)

Compute the spherical harmonic `Y_l^m(θ, φ)` numerically for scalar arguments.

The convention is

Y_l^m(θ, φ) = N_lm * P_l^m(cos(θ)) * exp(i m φ),

with

N_lm = sqrt((2l + 1)/(4π) * (l - |m|)! / (l + |m|)!).
"""
function spherical_harmonic(
    l::Integer,
    m::Integer,
    θ::Number,
    φ::Number,
    factorials::Union{Nothing,AbstractVector{<:Integer}}=nothing,
)::Number

    l < 0 && throw(ArgumentError("l must be non-negative"))
    abs(m) > l && throw(ArgumentError("|m| must be <= l"))
    facts = resolve_factorials(factorials, max(2l, abs(m) + l))
    N_lm = N(l, m, facts)
    P_lm = associated_legendre(l, m, cos(θ), factorials)
    Y_lm = N_lm * P(m, φ) * P_lm

    return Y_lm
end


"""
    spherical_harmonic(l, m, θ, φ)

Compute the spherical harmonic `Y_l^m(θ, φ)` symbolically for scalar arguments.

The convention is

Y_l^m(θ, φ) = N_lm * P_l^m(cos(θ)) * exp(i m φ),

with

N_lm = sqrt((2l + 1)/(4π) * (l - |m|)! / (l + |m|)!).
"""
function spherical_harmonic(
    l::Integer,
    m::Integer,
    θ::Num,
    φ::Num,
    factorials::Union{Nothing,AbstractVector{<:Integer}}=nothing,
)::Num

    l < 0 && throw(ArgumentError("l must be non-negative"))
    abs(m) > l && throw(ArgumentError("|m| must be <= l"))

    facts = resolve_factorials(factorials, max(2l, abs(m) + l))

    N_lm = N(l, m, facts)
    P_lm = associated_legendre(l, m, cos(θ), facts)
    Y_lm = N_lm * P(m, φ) * P_lm

    return expand(simplify(Y_lm))
end


"""
    spherical_harmonic(l, m, θ, φ; factorials=nothing)

Compute the spherical harmonic `Y_l^m(θ, φ)` elementwise on array-valued grids.
"""
function spherical_harmonic(
    l::Integer,
    m::Integer,
    θ::AbstractVector{<:Number},
    φ::AbstractVector{<:Number},
    factorials::Union{Nothing,AbstractVector{<:Integer}}=nothing,
)::AbstractMatrix{<:Number}
    facts = resolve_factorials(factorials, max(2l, abs(m) + l))
    N_lm = N(l, m, facts)
    P_lm = associated_legendre(l, m, cos.(θ), facts)   # length(θ)
    phase = P(m, φ)                                      # length(φ)

    P_lm_col = reshape(P_lm, :, 1)     # (length(θ), 1)
    phase_row = reshape(phase, 1, :)    # (1, length(φ))
    return N_lm .* P_lm_col .* phase_row   # -> (length(θ), length(φ)) macierz, bez wymogu równej długości
end
