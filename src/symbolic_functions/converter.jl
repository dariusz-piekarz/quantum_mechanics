using Symbolics

function _evaluate_numeric(expr, values)
    substituted = Symbolics.substitute(expr, values)
    numeric_expr = Symbolics.toexpr(substituted)
    return Base.eval(@__MODULE__, numeric_expr)
end

function to_function(expression, variables)
    vars = Tuple(collect(variables))

    scalar_fn = let expr = expression, vars = vars
        function scalar(args...)
            values = Dict(zip(vars, args))
            return _evaluate_numeric(expr, values)
        end
    end

    in_place_fn = let expr = expression, vars = vars
        function scalar!(out, args...)
            values = Dict(zip(vars, args))
            out[] = _evaluate_numeric(expr, values)
            return out
        end
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