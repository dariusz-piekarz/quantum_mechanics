using Symbolics

# ============================================================
# Legendre polynomial P_l(x)
# ============================================================

"""
    legendre(l, x)

Compute the Legendre polynomial `P_l(x)` numerically.
"""
function legendre(l::Integer, x::Number)::Number
    l < 0 && throw(ArgumentError("l must be non-negative"))

    l == 0 && return one(eltype(x))
    l == 1 && return x

    P₀ = one(eltype(x))
    P₁ = copy(x)

    for n in 2:l
        P₂ = ((2n - 1) * x * P₁ - (n - 1) * P₀) / n
        P₀, P₁ = P₁, P₂
    end

    return P₁
end


"""
    legendre(l, x)

Compute the Legendre polynomial `P_l(x)` symbolically.
"""
function legendre(l::Integer, x::Num)::Num
    l < 0 && throw(ArgumentError("l must be non-negative"))

    l == 0 && return one(eltype(x))
    l == 1 && return x

    P₀ = one(eltype(x))
    P₁ = copy(x)

    for n in 2:l
        P₂ = ((2n - 1) * x * P₁ - (n - 1) * P₀) / n
        P₀, P₁ = P₁, P₂
    end

    return simplify(expand(P₁))
end


"""
    legendre(l, x)

Evaluate the Legendre polynomial `P_l(x)` recursively on an array input.
"""
function legendre(l::Integer, x::AbstractArray{<:Number})::AbstractArray{<:Number}
    l < 0 && throw(ArgumentError("l must be non-negative"))

    # The recurrence divides by an integer, so integer arrays need a floating
    # output buffer to avoid InexactError on non-integral polynomial values.
    x_float = float.(x)
    l == 0 && return ones(eltype(x_float), size(x_float))
    l == 1 && return x_float

    P₀ = ones(eltype(x_float), size(x_float))
    P₁ = copy(x_float)
    buf = similar(x_float)

    for n in 2:l
        @. buf = ((2n - 1) * x_float * P₁ - (n - 1) * P₀) / n
        P₀, P₁, buf = P₁, buf, P₀
    end

    return P₁
end


# ============================================================
# Associated Legendre polynomial P_l^m(x)
# ============================================================

"""
    associated_legendre(l, m, x; factorials=nothing)

Compute the associated Legendre polynomial `P_l^m(x)` recursively over array data.
"""
function associated_legendre(
    l::Integer,
    m::Integer,
    x::AbstractArray{<:Number},
    factorials::Union{Nothing,AbstractVector{<:Integer}}=nothing,
)::AbstractArray{<:Number}

    l < 0 && throw(ArgumentError("l must be non-negative"))
    abs(m) > l && throw(ArgumentError("|m| must be <= l"))

    x_float = float.(x)
    ma = abs(m)

    if ma == 0
        Pmm = ones(eltype(x_float), size(x_float))
    else
        s = @. sqrt(1 - x_float^2)
        Pmm = ones(eltype(x_float), size(x_float))

        for k in 1:ma
            Pmm = @. -(2k - 1) * s * Pmm
        end
    end

    if l == ma
        result = Pmm
    elseif l == ma + 1
        result = @. (2ma + 1) * x_float * Pmm
    else
        Pprev = Pmm
        Pcurr = @. (2ma + 1) * x_float * Pmm
        buf = similar(x_float)

        for j in (ma+2):l
            @. buf = (
                (2j - 1) * x_float * Pcurr -
                (j + ma - 1) * Pprev
            ) / (j - ma)

            Pprev, Pcurr, buf = Pcurr, buf, Pprev
        end

        result = Pcurr
    end

    if m < 0
        facts = resolve_factorials(factorials, l + ma)
        coefficient = eltype(x_float)((-1)^ma * facts[l-ma+1] / facts[l+ma+1])
        result = @. coefficient * result
    end

    return result
end


"""
    associated_legendre(l, m, x)

Compute the associated Legendre polynomial `P_l^m(x)` symbolically.
"""
function associated_legendre(
    l::Integer,
    m::Integer,
    x::Num,
    factorials::Union{Nothing,AbstractVector{<:Integer}}=nothing,
)::Num

    l < 0 && throw(ArgumentError("l must be non-negative"))
    abs(m) > l && throw(ArgumentError("|m| must be <= l"))

    ma = abs(m)

    if ma == 0
        Pmm = one(eltype(x))
    else
        s = sqrt(1 - x^2)
        Pmm = one(eltype(x))

        for k in 1:ma
            Pmm = -(2k - 1) * s * Pmm
        end
    end

    if l == ma
        result = Pmm
    elseif l == ma + 1
        result = (2ma + 1) * x * Pmm
    else
        Pprev = Pmm
        Pcurr = (2ma + 1) * x * Pmm

        for j in (ma+2):l
            Pnext = (
                (2j - 1) * x * Pcurr -
                (j + ma - 1) * Pprev
            ) / (j - ma)

            Pprev, Pcurr = Pcurr, Pnext
        end

        result = Pcurr
    end

    if m < 0
        facts = resolve_factorials(factorials, l + ma)
        coefficient = (-1)^ma * facts[l-ma+1] / facts[l+ma+1]
        result = coefficient * result
    end

    return simplify(expand(result))
end


# Fallback version for scalar Number arguments
function associated_legendre(
    l::Integer,
    m::Integer,
    x::Number,
    factorials::Union{Nothing,AbstractVector{<:Integer}}=nothing,
)::Number

    l < 0 && throw(ArgumentError("l must be non-negative"))
    abs(m) > l && throw(ArgumentError("|m| must be <= l"))

    ma = abs(m)

    if ma == 0
        Pmm = one(eltype(x))
    else
        s = sqrt(1 - x^2)
        Pmm = one(eltype(x))

        for k in 1:ma
            Pmm = -(2k - 1) * s * Pmm
        end
    end

    if l == ma
        result = Pmm
    elseif l == ma + 1
        result = (2ma + 1) * x * Pmm
    else
        Pprev = Pmm
        Pcurr = (2ma + 1) * x * Pmm

        for j in (ma+2):l
            Pnext = (
                (2j - 1) * x * Pcurr -
                (j + ma - 1) * Pprev
            ) / (j - ma)

            Pprev, Pcurr = Pcurr, Pnext
        end

        result = Pcurr
    end

    if m < 0
        facts = resolve_factorials(factorials, l + ma)
        coefficient = (-1)^ma * facts[l-ma+1] / facts[l+ma+1]
        result = coefficient * result
    end

    return result
end
