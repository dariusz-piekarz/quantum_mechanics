using Symbolics

function to_function(expression, variables)
    variables = collect(variables)
    f = build_function(expression, variables...; expression=Val(false))
    if f isa Tuple
        return f[1], f[2]   # (out-of-place, in-place!)
    end
    return f, nothing
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