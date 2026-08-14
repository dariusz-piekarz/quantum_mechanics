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


# ============================================================
# Fibonacci numbers
# ============================================================


"""
Return the n-th Fibonacci number as a BigInt.

Arguments:
- n: non-negative integer index in the Fibonacci sequence

Returns:
- the n-th Fibonacci number using 1-based indexing (fibonacci(0) == 1)
"""
function fibonacci(n::Int)::BigInt
    if n == 0 || n == 1
        return 1
    end
    arr = BigInt[1, 1]
    for i = 2:n
        push!(arr, arr[i] + arr[i-1])
    end

    return arr[end]
end
