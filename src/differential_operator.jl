using Symbolics


struct DifferentialOperator{T}
    variable::T
end


D(var) = DifferentialOperator(var)


function (D::DifferentialOperator)(f)
    return expand_derivatives(Differential(D.variable)(f))
end


function gradient(f, variables)
    return [D(var)(f) for var in variables]
end


function gradient(f::SymbolicFunction, variables=f.variables)
    return SymbolicFunction(gradient(f.expression, variables), variables)
end


function divergence(f, variables)
    @assert length(f) == length(variables)
    return sum(D(variables[i])(f[i]) for i in eachindex(variables))
end


function divergence(f::SymbolicFunction, variables=f.variables)
    return SymbolicFunction(divergence(f.expression, variables), variables)
end


function laplacian(f,  variables)
    return divergence(gradient(f, variables), variables)
end


function laplacian(f::SymbolicFunction, variables=f.variables)
    return SymbolicFunction(laplacian(f.expression, variables),variables)
end


function hessian(f, variables)
    n = length(variables)
    h = Matrix{Any}(undef, n, n)

    for i in 1:n
        for j in 1:i
            h[i, j] = D(variables[i])(D(variables[j])(f))
            h[j, i] = h[i, j]
        end
    end
    
    return h
end


function hessian(f::SymbolicFunction, variables=f.variables)
    return SymbolicFunction(hessian(f.expression, variables), variables)
end




