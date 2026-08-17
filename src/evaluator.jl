using QuantumMechanics
using Symbolics


# ============================================================
# Scalar evaluation
# ============================================================

function _evaluate_1d!(f, X::AbstractMatrix, out::AbstractVector)
    @inbounds @simd for i in axes(X, 1)
        out[i] = f(X[i, 1])
    end

    return out
end


function _evaluate_2d!(f, X::AbstractMatrix, out::AbstractVector)
    @inbounds @simd for i in axes(X, 1)
        out[i] = f(
            X[i, 1],
            X[i, 2]
        )
    end

    return out
end


function _evaluate_3d!(f, X::AbstractMatrix, out::AbstractVector)
    @inbounds @simd for i in axes(X, 1)
        out[i] = f(
            X[i, 1],
            X[i, 2],
            X[i, 3]
        )
    end

    return out
end


function _evaluate_4d!(f, X::AbstractMatrix, out::AbstractVector)
    @inbounds @simd for i in axes(X, 1)
        out[i] = f(
            X[i, 1],
            X[i, 2],
            X[i, 3],
            X[i, 4]
        )
    end

    return out
end


function _evaluate_scalar_generic!(f, X::AbstractMatrix, out::AbstractVector)
    @inbounds @simd for i in axes(X, 1)
        out[i] = f(X[i, :]...)
    end

    return out
end


# ============================================================
# Array evaluation
#
# Result layout:
#
# scalar:
#     (npoints)
#
# vector:
#     (ncomponents, npoints)
#
# matrix:
#     (n1, n2, npoints)
#
# etc.
# ============================================================

function _evaluate_array_1d!(f, X::AbstractMatrix, out::AbstractArray)
    out_indices = ntuple(_ -> Colon(), ndims(out) - 1)
    @inbounds @simd for i in axes(X, 1)
        result = f(X[i, 1])

        out[out_indices..., i] .= result
    end

    return out
end


function _evaluate_array_2d!(f, X::AbstractMatrix, out::AbstractArray)
    out_indices = ntuple(_ -> Colon(), ndims(out) - 1)
    @inbounds @simd for i in axes(X, 1)
        result = f(
            X[i, 1],
            X[i, 2]
        )

        out[out_indices..., i] .= result
    end

    return out
end


function _evaluate_array_3d!(f, X::AbstractMatrix, out::AbstractArray)
    out_indices = ntuple(_ -> Colon(), ndims(out) - 1)
    @inbounds @simd for i in axes(X, 1)
        result = f(
            X[i, 1],
            X[i, 2],
            X[i, 3]
        )

        out[out_indices..., i] .= result
    end

    return out
end


function _evaluate_array_4d!(f, X::AbstractMatrix, out::AbstractArray)
    out_indices = ntuple(_ -> Colon(), ndims(out) - 1)

    @inbounds @simd for i in axes(X, 1)
        result = f(
            X[i, 1],
            X[i, 2],
            X[i, 3],
            X[i, 4]
        )

        out[out_indices..., i] .= result
    end

    return out
end


function _evaluate_array_generic!(f,X::AbstractMatrix,out::AbstractArray)
    npoints = size(X, 1)

    out_indices = ntuple(
        _ -> Colon(),
        ndims(out) - 1
    )

    @inbounds @simd for i in 1:npoints
        result = f(X[i, :]...)

        out[out_indices..., i] .= result
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

    # ========================================================
    # Scalar function
    # ========================================================

    if ndims(out) == 1

        fn = f.function_object

        if nvariables == 1

            return _evaluate_1d!(fn, X, out)

        elseif nvariables == 2

            return _evaluate_2d!(fn, X, out)

        elseif nvariables == 3

            return _evaluate_3d!(fn, X, out)

        elseif nvariables == 4

            return _evaluate_4d!(fn, X, out)

        else
            return _evaluate_scalar_generic!(fn, X, out)
        end

    # ========================================================
    # Non-scalar function
    # ========================================================

    else

        fn = f.function_object

        if nvariables == 1

            return _evaluate_array_1d!(fn, X, out)

        elseif nvariables == 2

            return _evaluate_array_2d!(fn, X, out)

        elseif nvariables == 3

            return _evaluate_array_3d!(fn,  X, out)

        elseif nvariables == 4

            return _evaluate_array_4d!(fn, X, out)

        else

            return _evaluate_array_generic!(fn, X, out)
        end
    end
end


# ============================================================
# evaluate
# ============================================================

function evaluate(f::SymbolicFunction, X::AbstractMatrix, ::Type{T} = Float64) where {T}
    npoints = size(X, 1)

    out = Array{T}(undef, size(f.expression)..., npoints)

    return evaluate!(f, X, out)
end
