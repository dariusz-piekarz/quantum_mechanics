using Symbolics

include(joinpath(@__DIR__, "math_special.jl"))


# ============================================================
# Legendre polynomial P_l(x)
# ============================================================

"""
    legendre(l, x::Union{Num, Number})

Compute the Legendre polynomial P_l(x) symbolically or numerically.
"""
function legendre(l::Integer, x::Union{Num,Number})::Union{Num,Number}
    l < 0 && throw(ArgumentError("l must be non-negative"))

    fact = factorials(2l) # Precompute factorials for efficiency
    lower_bound = ceil(Int, l / 2)

    polynomials = [(-1)^(l - k) * fact[2k+1] / (fact[k+1] * fact[2k-l+1] * fact[l-k+1]) * x^(2k - l) for k in lower_bound:l]
    if x isa Num
        return simplify(expand(sum(polynomials) / 2^l))
    else
        return sum(polynomials) / 2^l
    end
end


# ============================================================
# Associated Legendre polynomial P_l^m(x)
# ============================================================

"""
    associated_legendre(l, m, x::Union{Num,Number})

Compute the associated Legendre polynomial P_l^m(x)
symbolically or numerically.
"""
function associated_legendre(l::Integer, m::Integer, x::Union{Num,Number})::Union{Num,Number}
    l < 0 && throw(ArgumentError("l must be non-negative"))
    abs(m) > l && throw(ArgumentError("|m| must be <= l"))
    m_abs = abs(m)
    lower_bound = ceil(Int, (m_abs + l) / 2)
    fact = factorials(2l) # Precompute factorials for efficiency
    polynomials = [
        (-1)^(l-k+m_abs)*fact[2k+1] / (fact[k+1] * fact[l-k+1] * fact[2k-l-m_abs+1]) * x^(2k-l-m_abs) for k in lower_bound:l
    ]
    result = sum(polynomials) * (1-x^2)^(m_abs / 2) / 2^l
    if m > 0 && x isa Num
        return simplify(expand(result))
    elseif m > 0
        return result
    elseif m < 0 && x isa Num
        return simplify(expand((-1)^m_abs * factorial(l-m_abs) / factorial(l+m_abs) * result))
    else
        return (-1)^m_abs * factorial(l-m_abs) / factorial(l+m_abs) * result
    end
end
