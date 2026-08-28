# ============================================================
# factorials
# ============================================================

"""
    factorials_table(n)

Return the factorial table for `0, 1, ..., n` as a vector, with the convention
`facts[i+1] = i!` (so `facts[1] = 0! = 1`).

`Int64` overflows already at `21!` (`21! ≈ 5.11e19 > typemax(Int64)`), which
would otherwise silently wrap around and corrupt every downstream calculation
(radial functions, associated Legendre polynomials, spherical harmonics) for
`n + l > 20`. To stay correct for higher quantum numbers while keeping the
common case fast, this returns `Vector{Int}` when it's safe (`n <= 20`) and
transparently switches to `Vector{BigInt}` otherwise.

# Arguments
- `n::Integer`: non-negative integer up to which factorials are computed.

# Returns
A vector of factorial values, `Vector{Int}` or `Vector{BigInt}` depending on `n`.
"""
function factorials_table(n::Integer)::Vector{<:Integer}
    n < 0 && throw(ArgumentError("n must be non-negative"))
    T = n <= 20 ? Int : BigInt
    facts = Vector{T}(undef, n + 1)
    facts[1] = one(T)
    for i in 2:(n+1)
        facts[i] = facts[i-1] * T(i - 1)
    end
    return facts
end

"""
    resolve_factorials(facts, n)

Return a factorial table covering `0, 1, ..., n`.
If `facts` is `nothing`, a new table is generated. Otherwise, the provided vector
is validated and reused. The expected convention is `facts[i+1] = i!`.
"""
function resolve_factorials(
    facts::Union{Nothing,AbstractVector{<:Integer}},
    n::Integer,
)
    if facts === nothing
        return factorials_table(n)
    end
    if length(facts) <= n
        throw(ArgumentError("factorials_table vector is too short for n = $n"))
    end
    return facts
end