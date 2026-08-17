using Symbolics


function to_function(expression, variables)

    variables = collect(variables)

    f = build_function(
        expression,
        variables...;
        expression = Val(false)
    )

    if f isa Tuple
        return first(f)
    end

    return f
end


struct SymbolicFunction{E,V,F}
    expression::E
    variables::V
    function_object::F
end


function SymbolicFunction(expression, variables)
    f = to_function(expression, variables)
    return SymbolicFunction(expression, variables, f)
end


(f::SymbolicFunction)(args...) = f.function_object(args...)
