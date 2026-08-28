using Symbolics
import ForwardDiff


_is_symbolic(variables) = eltype(variables) <: Num


struct DifferentialOperator{T}
    variable::T
end


D(var) = DifferentialOperator(var)


function (D::DifferentialOperator)(f)
    return expand_derivatives(Differential(D.variable)(f))
end

# ============================================================
# Gradient
# ============================================================

function gradient(f, variables)
    if _is_symbolic(variables)
        return [D(var)(f) for var in variables]
    end
    return ForwardDiff.gradient(f, variables)
end

# SymbolicFunction
function gradient(f::SymbolicFunction, variables=f.variables)
    return SymbolicFunction(gradient(f.expression, variables), variables)
end

# ============================================================
# Divergence
# ============================================================

function divergence(f, variables)
    
    if _is_symbolic(variables)
        @assert length(f) == length(variables)
        
        return sum(D(variables[i])(f[i]) for i in eachindex(variables))
    end

    jacobian = ForwardDiff.jacobian(f, variables)
    return sum(jacobian[i, i] for i in 1: length(jacobian[:, 1]))
end


# SymbolicFunction
function divergence(f::SymbolicFunction,variables=f.variables)
    return SymbolicFunction(divergence(f.expression, variables), variables)
end

# ============================================================
# Laplacian
# ============================================================

function laplacian(f, variables)
    if _is_symbolic(variables)
        return divergence(gradient(f, variables),variables)
    end

    hessian = ForwardDiff.hessian(f, variables)

    return sum(hessian[i, i] for i in 1: length(hessian[:, 1]))
end

# SymbolicFunction
function laplacian(f::SymbolicFunction, variables=f.variables)
    return SymbolicFunction(laplacian(f.expression, variables), variables)
end

# ============================================================
# Hessian
# ============================================================

function hessian(f, variables)
    if _is_symbolic(variables)
         n = length(variables)

        h = Matrix{Num}(undef, n, n)

        for i in 1:n
            for j in 1:i
                h[i, j] =
                    D(variables[i])(
                        D(variables[j])(f)
                    )

                h[j, i] = h[i, j]
            end
        end

        return h
    end

    return ForwardDiff.hessian(f, variables)
end

# SymbolicFunction
function hessian(f::SymbolicFunction, variables=f.variables)
    return SymbolicFunction(hessian(f.expression, variables), variables)
end
