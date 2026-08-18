using QuantumMechanics

abstract type AbstractOperator end


struct KineticOperator{T} <: AbstractOperator
    coefficient::T
end


function KineticOperator(::Type{T}) where {T}
    return KineticOperator{T}(-one(T) / T(2))
end


function (T̂::KineticOperator{T})(f::SymbolicFunction) where {T}
    lap = laplacian(f)

    expression = T̂.coefficient * lap.expression

    return SymbolicFunction(expression, f.variables)
end