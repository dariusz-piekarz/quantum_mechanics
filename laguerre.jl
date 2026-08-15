using Symbolics


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
function laguerre(n::Integer, x::Number, factorials::Union{Nothing,AbstractVector{<:Integer}}=nothing)::Number
    n < 0 && throw(ArgumentError("n must be non-negative"))

    fact = resolve_factorials(factorials, n)
    polynomials = [(-1)^k * fact[n+1] / (fact[k+1]^2 * fact[n-k+1]) * x^k for k in 0:n]
    return sum(polynomials)

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
function laguerre(n::Integer, x::Num, factorials::Union{Nothing,AbstractVector{<:Integer}}=nothing)::Num
    n < 0 && throw(ArgumentError("n must be non-negative"))

    fact = resolve_factorials(factorials, n)
    polynomials = [(-1)^k * fact[n+1] / (fact[k+1]^2 * fact[n-k+1]) * x^k for k in 0:n]

    return simplify(expand(sum(polynomials)))

end


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
function generalized_laguerre(
    n::Integer,
    α::Number,
    x::Number;
    factorials::Union{Nothing,AbstractVector{<:Integer}}=nothing
)::Number

    n < 0 && throw(ArgumentError("n must be non-negative"))
    if !(α isa Num) && α < -1
        throw(ArgumentError("α must be >= -1"))
    end
    # L_n^α(x) = Σ(k=0 to n) [(-1)^k * Γ(n+α+1) / (Γ(α+k+1) * Γ(n-k+1) * k!)] * x^k
    # Simplified: = Σ(k=0 to n) [(-1)^k / k! * prod((α+k+1):(α+n))] * x^k
    fact = resolve_factorials(factorials, n)
    polynomials = []

    for k in 0:n
        # Compute the rising product (α+k+1) * (α+k+2) * ... * (α+n)
        prod_val = prod((α + j) for j in (k+1):n; init=1)
        coeff = (-1)^k * prod_val / (fact[k+1] * fact[n-k+1])
        push!(polynomials, coeff * x^k)
    end

    result = sum(polynomials)
    return result
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
function generalized_laguerre(
    n::Integer,
    α::Union{Num,Number},
    x::Num;
    factorials::Union{Nothing,AbstractVector{<:Integer}}=nothing
)::Num

    n < 0 && throw(ArgumentError("n must be non-negative"))
    if ~(α isa Num) && α < -1
        throw(ArgumentError("α must be >= -1"))
    end
    # L_n^α(x) = Σ(k=0 to n) [(-1)^k * Γ(n+α+1) / (Γ(α+k+1) * Γ(n-k+1) * k!)] * x^k
    # Simplified: = Σ(k=0 to n) [(-1)^k / k! * prod((α+k+1):(α+n))] * x^k
    fact = resolve_factorials(factorials, n)

    polynomials = []

    for k in 0:n
        # Compute the rising product (α+k+1) * (α+k+2) * ... * (α+n)
        prod_val = 1
        for j in (k+1):n
            prod_val = prod_val * (α + j)
        end
        coeff = (-1)^k * prod_val / (fact[k+1] * fact[n-k+1])
        push!(polynomials, coeff * x^k)
    end

    result = sum(polynomials)

    return simplify(expand(result))
end

# Backward-compatible alias
function generalized_lagguerre(
    n::Integer,
    α,
    x;
    factorials::Union{Nothing,AbstractVector{<:Integer}}=nothing,
)
    return generalized_laguerre(n, α, x; factorials=factorials)
end


"""
    generalized_laguerre(n, α, x)

Compute the generalized Laguerre polynomial `L_n^(α)(x)` on array inputs.
"""
function generalized_laguerre(
    n::Integer,
    α::Number,
    x::AbstractArray{<:Number}=nothing
)
    n < 0 && throw(ArgumentError("n must be non-negative"))
    !(α isa Num) && α < -1 && throw(ArgumentError("α must be >= -1"))
    x === nothing && throw(ArgumentError("x must not be nothing"))

    L0 = ones(eltype(x), size(x))

    n == 0 && return L0

    L1 = @. 1 + α - x

    n == 1 && return L1

    L_prev = L0
    L_curr = L1

    for k in 2:n
        L_next = @. (
            (2k - 1 + α - x) * L_curr -
            (k - 1 + α) * L_prev
        ) / k

        L_prev, L_curr = L_curr, L_next
    end

    return L_curr
end
