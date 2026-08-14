using Symbolics

include(joinpath(@__DIR__, "math_special.jl"))


"""
    laguerre(n::Integer, x::Union{Num,Number})::Union{Num,Number}

Compute the Laguerre polynomial L_n(x) of degree n.

The Laguerre polynomials are orthogonal polynomials used in quantum mechanics,
particularly in solutions of the hydrogen atom.

# Arguments
- `n::Integer`: Degree of the polynomial (must be non-negative)
- `x::Union{Num,Number}`: Evaluation point (can be symbolic or numeric)

# Returns
- `Union{Num,Number}`: The value of the Laguerre polynomial L_n(x)

# Throws
- `ArgumentError`: If n < 0

# Examples
```julia
using Symbolics
@syms x
laguerre(2, x)  # Returns the 2nd degree Laguerre polynomial
laguerre(3, 1.5)  # Returns numerical value
```
"""
function laguerre(n::Integer, x::Union{Num,Number})::Union{Num,Number}
    n < 0 && throw(ArgumentError("n must be non-negative"))

    fact = factorials(n) # Precompute factorials for efficiency
    polynomials = [(-1)^k * fact[n+1] / (fact[k+1]^2 * fact[n-k+1]) * x^k for k in 0:n]
    if x isa Num
        return simplify(expand(sum(polynomials)))
    else
        return sum(polynomials)
    end
end


"""
    generalized_lagguerre(n::Integer, α::Union{Num,Number}, x::Union{Num,Number})::Union{Num,Number}

Compute the generalized Laguerre polynomial L_n^(α)(x) of degree n with parameter α.

The generalized Laguerre polynomials are a family of orthogonal polynomials parametrized
by α. They reduce to standard Laguerre polynomials when α = 0.

# Arguments
- `n::Integer`: Degree of the polynomial (must be non-negative)
- `α::Union{Num,Number}`: Parameter for generalization
- `x::Union{Num,Number}`: Evaluation point (can be symbolic or numeric)

# Returns
- `Union{Num,Number}`: The value of the generalized Laguerre polynomial L_n^(α)(x)

# Throws
- `ArgumentError`: If n < 0 or α < -1

# Examples
```julia
using Symbolics
@syms x
generalized_lagguerre(2, 0, x)  # Standard Laguerre polynomial
generalized_lagguerre(1, 0.5, 2.0)  # Numerical evaluation with α = 0.5
```
"""
function generalized_lagguerre(n::Integer, α::Union{Num,Number}, x::Union{Num,Number})::Union{Num,Number}
    n < 0 && throw(ArgumentError("n must be non-negative"))
    if ~(α isa Num) && α < -1
        throw(ArgumentError("α must be >= -1"))
    end
    # L_n^α(x) = Σ(k=0 to n) [(-1)^k * Γ(n+α+1) / (Γ(α+k+1) * Γ(n-k+1) * k!)] * x^k
    # Simplified: = Σ(k=0 to n) [(-1)^k / k! * prod((α+k+1):(α+n))] * x^k

    fact = factorials(n) # Precompute factorials for efficiency
    polynomials = []

    for k in 0:n
        # Compute the rising product (α+k+1) * (α+k+2) * ... * (α+n)
        if ~(α isa Num)
            prod_val = prod((α+k+1):(α+n))
        else
            # For symbolic α, compute the product symbolically
            prod_val = 1
            for j in (k+1):n
                prod_val = prod_val * (α + j)
            end
        end

        coeff = (-1)^k * prod_val / (fact[k+1] * fact[n-k+1])
        push!(polynomials, coeff * x^k)
    end

    result = sum(polynomials)
    if x isa Num
        return simplify(expand(result))
    else
        return result
    end

end
