using Symbolics

# ============================================================
# Laguerre polynomial L_n(x)
# ============================================================

"""
    laguerre(n::Integer, x::Number)::Number

Compute the Laguerre polynomial L_n(x) of degree n.

The Laguerre polynomials are orthogonal polynomials used in quantum mechanics,
particularly in solutions of the hydrogen atom.

# Arguments
- `n::Integer`: Degree of the polynomial (must be non-negative)
- `x::Number`: Evaluation point (numeric)

# Returns
- `Number`: The value of the Laguerre polynomial L_n(x)

# Throws
- `ArgumentError`: If n < 0

# Examples
```julia
laguerre(3, 1.5)  # Returns numerical value
```
"""
function laguerre(n::Integer, x::Number)::Number
    n < 0 && throw(ArgumentError("n must be non-negative"))

    L0 = one(eltype(x))

    n == 0 && return L0

    L1 = 1 - x

    n == 1 && return L1

    L_prev = L0
    L_curr = L1

    for k in 2:n
        L_next = ((2k - 1 - x) * L_curr - (k - 1) * L_prev) / k
        L_prev, L_curr = L_curr, L_next
    end

    return L_curr

end


"""
    laguerre(n::Integer, x::Num)::Num

Compute the Laguerre polynomial L_n(x) of degree n.

The Laguerre polynomials are orthogonal polynomials used in quantum mechanics,
particularly in solutions of the hydrogen atom.

# Arguments
- `n::Integer`: Degree of the polynomial (must be non-negative)
- `x::Num`: Evaluation point (symbolic)

# Returns
- `Num`: The value of the Laguerre polynomial L_n(x)

# Throws
- `ArgumentError`: If n < 0

# Examples
```julia
using Symbolics
@syms x
laguerre(2, x)  # Returns the 2nd degree Laguerre polynomial
```
"""
function laguerre(n::Integer, x::Num)::Num
    n < 0 && throw(ArgumentError("n must be non-negative"))

    L0 = one(eltype(x))

    n == 0 && return L0

    L1 = 1 - x

    n == 1 && return L1

    L_prev = L0
    L_curr = L1

    for k in 2:n
        L_next = ((2k - 1 - x) * L_curr - (k - 1) * L_prev) / k
        L_prev, L_curr = L_curr, L_next
    end

    return simplify(expand(result))

end


# Backward-compatible alias
function lagguerre(n::Integer, x)
    return generalized_laguerre(n, x)
end


# ============================================================
# Generalized Laguerre polynomial L_n^(α)(x)
# ============================================================

"""
    generalized_laguerre(n::Integer, α::Number, x::Number)::Number

Compute the generalized Laguerre polynomial L_n^(α)(x) of degree n with parameter α.

The generalized Laguerre polynomials are a family of orthogonal polynomials parametrized
by α. They reduce to standard Laguerre polynomials when α = 0.

# Arguments
- `n::Integer`: Degree of the polynomial (must be non-negative)
- `α::Number`: Parameter for generalization
- `x::Number`: Evaluation point (numeric)

# Returns
- `Number`: The value of the generalized Laguerre polynomial L_n^(α)(x)

# Throws
- `ArgumentError`: If n < 0 or α < -1

# Examples
```julia
generalized_laguerre(1, 0.5, 2.0)  # Numerical evaluation with α = 0.5
```
"""
function generalized_laguerre(n::Integer, α::Number, x::Number)::Number
    n < 0 && throw(ArgumentError("n must be non-negative"))
    α < -1 && throw(ArgumentError("α must be >= -1"))

    L0 = one(eltype(x))

    n == 0 && return L0

    L1 = 1 + α - x

    n == 1 && return L1

    L_prev = L0
    L_curr = L1

    for k in 2:n
        L_next = ((2k - 1 + α - x) * L_curr - (k - 1 + α) * L_prev) / k
        L_prev, L_curr = L_curr, L_next
    end

    return L_curr
end


"""
    generalized_laguerre(n::Integer, α::Union{Num,Number}, x::Num)::Num

Compute the generalized Laguerre polynomial L_n^(α)(x) of degree n with parameter α.

The generalized Laguerre polynomials are a family of orthogonal polynomials parametrized
by α. They reduce to standard Laguerre polynomials when α = 0.

# Arguments
- `n::Integer`: Degree of the polynomial (must be non-negative)
- `α::Union{Num,Number}`: Parameter for generalization
- `x::Num`: Evaluation point (symbolic)

# Returns
- `Num`: The value of the generalized Laguerre polynomial L_n^(α)(x)

# Throws
- `ArgumentError`: If n < 0 or α < -1

# Examples
```julia
using Symbolics
@syms x
generalized_laguerre(2, 0, x)  # Standard Laguerre polynomial
```
"""
function generalized_laguerre(n::Integer, α::Union{Num,Number}, x::Num)::Num

    n < 0 && throw(ArgumentError("n must be non-negative"))
    if ~(α isa Num) && α < -1
        throw(ArgumentError("α must be >= -1"))
    end

    L0 = one(eltype(x))

    n == 0 && return L0

    L1 = 1 + α - x

    n == 1 && return L1

    L_prev = L0
    L_curr = L1

    for k in 2:n
        L_next = ((2k - 1 + α - x) * L_curr - (k - 1 + α) * L_prev) / k
        L_prev, L_curr = L_curr, L_next
    end

    return simplify(expand(L_curr))
end


"""
    generalized_laguerre(n, α, x)

Compute the generalized Laguerre polynomial `L_n^(α)(x)` on array inputs.
"""
function generalized_laguerre(
    n::Integer,
    α::Number,
    x::AbstractArray{<:Number}
)::AbstractArray{<:Number}
    n < 0 && throw(ArgumentError("n must be non-negative"))
    α < -1 && throw(ArgumentError("α must be >= -1"))

    L0 = ones(eltype(x), size(x))

    n == 0 && return L0

    L1 = @. 1 + α - x

    n == 1 && return L1

    L_prev = L0
    L_curr = L1

    for k in 2:n
        L_next = @. ((2k - 1 + α - x) * L_curr - (k - 1 + α) * L_prev) / k
        L_prev, L_curr = L_curr, L_next
    end

    return L_curr
end


# Backward-compatible alias
function generalized_lagguerre(n::Integer, α, x)
    return generalized_laguerre(n, α, x)
end

