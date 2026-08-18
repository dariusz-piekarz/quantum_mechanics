using Symbolics
import ForwardDiff


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

# Symbolic
function gradient(f, variables)
    return [D(var)(f) for var in variables]
end

# Numerical
function gradient(f, variables::AbstractVector)
    return ForwardDiff.gradient(f, variables)
end

# SymbolicFunction
function gradient(f::SymbolicFunction, variables=f.variables)
    return SymbolicFunction(gradient(f.expression, variables), variables)
end


# ============================================================
# Divergence
# ============================================================

# Symbolic
function divergence(f, variables)
    @assert length(f) == length(variables)

    return sum(
        D(variables[i])(f[i])
        for i in eachindex(variables)
    )
end

# Numerical
function divergence(f, variables::AbstractVector)
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

# Symbolic
function laplacian(f, variables)
    return divergence(gradient(f, variables),variables)
end

# Numerical
function laplacian(f, variables::AbstractVector)
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

# Symbolic
function hessian(f, variables)
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

# Numerical
function hessian(f, variables::AbstractVector)
    return ForwardDiff.hessian(f, variables)
end

# SymbolicFunction
function hessian(f::SymbolicFunction, variables=f.variables)
    return SymbolicFunction(hessian(f.expression, variables), variables)
end
