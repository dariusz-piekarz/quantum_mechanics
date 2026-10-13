using Symbolics

function _check_arity(args, nvariables)
    length(args) == nvariables ||
        throw(ArgumentError("expected $nvariables arguments, got $(length(args))"))
end

function to_function(expression, variables)
    vars = Tuple(collect(variables))
    nvariables = length(vars)

    if expression isa AbstractArray
        compiled, compiled! = Symbolics.build_function(
            expression, vars...; expression=Val(false), skipzeros=false
        )

        scalar_fn = function (args...)
            _check_arity(args, nvariables)
            return compiled(args...)
        end

        in_place_fn = function (out, args...)
            _check_arity(args, nvariables)
            if size(out) == size(expression)
                compiled!(out, args...)
            elseif length(out) == length(expression)
                compiled!(reshape(out, size(expression)), args...)
            else
                throw(DimensionMismatch("expected output size $(size(expression)), got $(size(out))"))
            end
            return out
        end

        return scalar_fn, in_place_fn
    end

    compiled = Symbolics.build_function(expression, vars...; expression=Val(false))
    scalar_fn = function (args...)
        _check_arity(args, nvariables)
        return compiled(args...)
    end
    in_place_fn = function (out, args...)
        _check_arity(args, nvariables)
        out .= compiled(args...)
        return out
    end
    return scalar_fn, in_place_fn
end

struct SymbolicFunction{E,V,F,FIP}
    expression::E
    variables::V
    function_object::F
    function_object!::FIP
end

function SymbolicFunction(expression, variables)
    f, fip = to_function(expression, variables)
    return SymbolicFunction(expression, variables, f, fip)
end

(f::SymbolicFunction)(args...) = f.function_object(args...)

function (f::SymbolicFunction)(X::AbstractMatrix)
    return evaluate(f, X)
end
