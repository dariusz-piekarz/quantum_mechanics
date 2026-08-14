include(joinpath(@__DIR__, "legendre.jl"))


# ============================================================
# Spherical harmonics Y_l^m(θ, φ)
# ============================================================

"""
    spherical_harmonic(l, m, θ::Num, φ::Num)

Compute the spherical harmonic Y_l^m(θ, φ) symbolically.

The convention is

Y_l^m(θ, φ)
=
N_lm P_l^m(cos(θ)) exp(i m φ)

where

N_lm =
sqrt(
    (2l + 1)/(4π)
    * (l-m)!/(l+m)!
)
"""
function spherical_harmonic(l::Integer, m::Integer, θ::Union{Num,Number}, φ::Union{Num,Number})

    l < 0 && throw(ArgumentError("l must be non-negative"))
    abs(m) > l && throw(ArgumentError("|m| must be <= l"))

    m_abs = abs(m)

    N_lm = sqrt((2l + 1) / (4π) * factorial(l - m_abs) / factorial(l + m_abs))

    P_lm = associated_legendre(l, m, cos(θ))

    Y_lm = N_lm * P_lm * exp(im * m * φ)

    if θ isa Num || φ isa Num
        return expand(simplify(Y_lm))
    else
        return Y_lm
    end
end


"""
    spherical_harmonic(l, m, θ, φ)

Compute Y_l^m elementwise for vectors θ and φ.
"""
function spherical_harmonic(
    l::Integer,
    m::Integer,
    θ::AbstractVector{<:Number},
    φ::AbstractVector{<:Number},
)

    return [spherical_harmonic(l, m, θᵢ, φⱼ) for θᵢ in θ, φⱼ in φ]
end
