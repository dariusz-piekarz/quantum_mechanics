using QuantumMechanics


function Hamiltonian(V::CoulombPotential{T}, K::KineticOperator{T}, ψ::SymbolicFunction) where {T}
    variables = ψ.variables

    Tψ = K(ψ)
    Vψ = V(variables...)

    expression = Tψ.expression + Vψ * ψ.expression

    return SymbolicFunction(expression, variables)
end