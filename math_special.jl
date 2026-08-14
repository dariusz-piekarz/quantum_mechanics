# ============================================================
# factorials
# ============================================================

"""
    factorials(n)

Return the factorial values for 0, 1, ..., n as a vector.

# Arguments
- `n::Integer`: non-negative integer up to which factorials are computed.

# Returns
A vector `facts` such that `facts[k+1] = k!` for `k = 0,1,...,n`.
"""
function factorials(n::Integer)::Vector{Int}
    n < 0 && throw(ArgumentError("n must be non-negative"))

    facts = Vector{Int}(undef, n + 1)
    facts[1] = 1

    for i in 2:(n+1)
        facts[i] = facts[i-1] * (i - 1)
    end

    return facts

end
