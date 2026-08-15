# ============================================================
# factorials
# ============================================================

"""
    factorials(n)

Return the factorial table for `0, 1, ..., n` as a vector.

The convention is `facts[i+1] = i!` for `i = 0, 1, ..., n`, so `facts[1] = 0! = 1`.

# Arguments
- `n::Integer`: non-negative integer up to which factorials are computed.

# Returns
A vector of integer factorial values.
"""
function factorials_table(n::Integer)::Vector{Int}
    n < 0 && throw(ArgumentError("n must be non-negative"))

    facts = Vector{Int}(undef, n + 1)
    facts[1] = 1

    for i in 2:(n+1)
        facts[i] = facts[i-1] * (i - 1)
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
