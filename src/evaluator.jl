using QuantumMechanics
using Symbolics

# ============================================================
# Scalar evaluation (f: out-of-place, out::AbstractVector)
# ============================================================

function _evaluate_1d!(f, X::AbstractMatrix, out::AbstractVector)
    @inbounds @simd for i in axes(X, 1)
        out[i] = f(X[i, 1])
    end
    return out
end

function _evaluate_2d!(f, X::AbstractMatrix, out::AbstractVector)
    @inbounds @simd for i in axes(X, 1)
        out[i] = f(X[i, 1], X[i, 2])
    end
    return out
end

function _evaluate_3d!(f, X::AbstractMatrix, out::AbstractVector)
    @inbounds @simd for i in axes(X, 1)
        out[i] = f(X[i, 1], X[i, 2], X[i, 3])
    end
    return out
end

function _evaluate_4d!(f, X::AbstractMatrix, out::AbstractVector)
    @inbounds @simd for i in axes(X, 1)
        out[i] = f(X[i, 1], X[i, 2], X[i, 3], X[i, 4])
    end
    return out
end

function _evaluate_scalar_generic!(f, X::AbstractMatrix, out::AbstractVector)
    @inbounds @simd for i in axes(X, 1)
        out[i] = f(view(X, i, :)...)
    end
    return out
end

# ============================================================
# Array evaluation (f!: in-place, out::AbstractArray)
#
# Result layout:
#   scalar expr -> out :: (npoints,)
#   vector expr -> out :: (ncomponents, npoints)
#   matrix expr -> out :: (n1, n2, npoints)
# ============================================================

function _evaluate_array_1d!(f!, X::AbstractMatrix, out::AbstractArray)
    out_indices = ntuple(_ -> Colon(), ndims(out) - 1)
    @inbounds @simd for i in axes(X, 1)
        f!(view(out, out_indices..., i), X[i, 1])
    end
    return out
end

function _evaluate_array_2d!(f!, X::AbstractMatrix, out::AbstractArray)
    out_indices = ntuple(_ -> Colon(), ndims(out) - 1)
    @inbounds @simd for i in axes(X, 1)
        f!(view(out, out_indices..., i), X[i, 1], X[i, 2])
    end
    return out
end

function _evaluate_array_3d!(f!, X::AbstractMatrix, out::AbstractArray)
    out_indices = ntuple(_ -> Colon(), ndims(out) - 1)
    @inbounds @simd for i in axes(X, 1)
        f!(view(out, out_indices..., i), X[i, 1], X[i, 2], X[i, 3])
    end
    return out
end

function _evaluate_array_4d!(f!, X::AbstractMatrix, out::AbstractArray)
    out_indices = ntuple(_ -> Colon(), ndims(out) - 1)
    @inbounds @simd for i in axes(X, 1)
        f!(view(out, out_indices..., i), X[i, 1], X[i, 2], X[i, 3], X[i, 4])
    end
    return out
end

function _evaluate_array_generic!(f!, X::AbstractMatrix, out::AbstractArray)
    out_indices = ntuple(_ -> Colon(), ndims(out) - 1)
    @inbounds @simd for i in axes(X, 1)
        f!(view(out, out_indices..., i), view(X, i, :)...)
    end
    return out
end

# ============================================================
# Main evaluate!
# ============================================================

function evaluate!(f::SymbolicFunction, X::AbstractMatrix, out::AbstractArray)
    npoints = size(X, 1)
    @assert size(out, ndims(out)) == npoints

    nvariables = length(f.variables)

    if ndims(out) == 1
        fn = f.function_object
        nvariables == 1 && return _evaluate_1d!(fn, X, out)
        nvariables == 2 && return _evaluate_2d!(fn, X, out)
        nvariables == 3 && return _evaluate_3d!(fn, X, out)
        nvariables == 4 && return _evaluate_4d!(fn, X, out)
        return _evaluate_scalar_generic!(fn, X, out)
    else
        fn! = f.function_object!
        nvariables == 1 && return _evaluate_array_1d!(fn!, X, out)
        nvariables == 2 && return _evaluate_array_2d!(fn!, X, out)
        nvariables == 3 && return _evaluate_array_3d!(fn!, X, out)
        nvariables == 4 && return _evaluate_array_4d!(fn!, X, out)
        return _evaluate_array_generic!(fn!, X, out)
    end
end

# ============================================================
# evaluate
# ============================================================

function evaluate(f::SymbolicFunction, X::AbstractMatrix{S}, ::Type{T}=real(S)) where {S,T}
    npoints = size(X, 1)
    out = Array{T}(undef, size(f.expression)..., npoints)
    return evaluate!(f, X, out)
end