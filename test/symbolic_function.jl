using Test
using Symbolics
using QuantumMechanics


@testset "SymbolicFunction" begin

    @variables x y z t

    # ---------------------------------------------------------
    # Construction
    # ---------------------------------------------------------

    @testset "construction" begin

        expr = x^2 + y^2 + z^2
        f = SymbolicFunction(expr, (x, y, z))

        @test isequal(f.expression, expr)
        @test f.variables == (x, y, z)
    end

    # ---------------------------------------------------------
    # Evaluation
    # ---------------------------------------------------------

    @testset "evaluation" begin

        expr = x^2 + y^2 + z^2
        f = SymbolicFunction(expr, (x, y, z))

        @test f(1.0, 2.0, 3.0) ≈ 14.0
        @test f(0.0, 0.0, 0.0) ≈ 0.0
        @test f(2.0, -3.0, 4.0) ≈ 29.0
    end

    # ---------------------------------------------------------
    # to_function
    # ---------------------------------------------------------

    @testset "to_function" begin

        expr = x^2 + y^2 + z^2

        f = to_function(expr, (x, y, z))

        @test f(1.0, 2.0, 3.0) ≈ 14.0
        @test f(2.0, -3.0, 4.0) ≈ 29.0
    end

    # ---------------------------------------------------------
    # Single variable
    # ---------------------------------------------------------

    @testset "single variable" begin
        expr = t^3 + 2t + 1

        f = SymbolicFunction(expr, (t,))

        @test f(0.0) ≈ 1.0
        @test f(1.0) ≈ 4.0
        @test f(2.0) ≈ 13.0

        g = to_function(expr, (t,))

        @test g(0.0) ≈ 1.0
        @test g(2.0) ≈ 13.0
    end


        # ---------------------------------------------------------
    # Mathematical functions
    # ---------------------------------------------------------

    @testset "elementary functions" begin

        expr = sin(x) + cos(y) + exp(z)

        f = SymbolicFunction(expr, (x, y, z))

        expected = sin(1.0) + cos(2.0) + exp(3.0)

        @test f(1.0, 2.0, 3.0) ≈ expected
    end

    @testset "numeric divergence/laplacian match symbolic" begin
    fvec = [x, y^2, z^3]
    expr = x^2 + y^2 + z^2

    div_sym = divergence(fvec, (x, y, z))     
    lap_sym = laplacian(expr, (x, y, z))

    point = [1.0, 2.0, 3.0]
    fvec_num(v) = [v[1], v[2]^2, v[3]^3]
    fscal_num(v) = v[1]^2 + v[2]^2 + v[3]^2

    @test divergence(fvec_num, point) ≈ 1 + 2*2.0 + 3*3.0^2
    @test laplacian(fscal_num, point) ≈ 6.0
end

end


@testset "Symbolic Differential Operators" begin
    @variables x y z
    # ---------------------------------------------------------
    # Gradient
    # ---------------------------------------------------------

    @testset "gradient" begin

        expr = x^2 + y^2 + z^2

        g = gradient(expr, (x,y,z))
        f = SymbolicFunction(g, (x, y, z))
        @test isequal((f.expression), [2x, 2y, 2z])

        result = f(1.0, 2.0, 3.0)
        @test result ≈ [2.0, 4.0, 6.0]
    end

    # ---------------------------------------------------------
    # Laplacian
    # ---------------------------------------------------------

    @testset "laplacian" begin

        expr = x^2 + y^2 + z^2

        lp = laplacian(expr, (x, y, z))
        Δf = SymbolicFunction(lp, (x, y, z))

        @test isequal(Δf.expression, 6)
        @test Δf(1.0, 2.0, 3.0) ≈ 6.0
        @test Δf(-10.0, 4.0, 7.0) ≈ 6.0
    end

    # ---------------------------------------------------------
    # Hessian
    # ---------------------------------------------------------

    @testset "hessian" begin

        expr = x^2 + y^2 + z^2

        f = hessian(expr, (x, y, z))
        H = SymbolicFunction(f, (x, y, z))

        expected = [
            2  0  0
            0  2  0
            0  0  2
        ]

        @test isequal(H.expression, expected)
        result = H(1.0, 2.0, 3.0)
        @test result ≈ Float64.(expected)
    end

        # ---------------------------------------------------------
    # Divergence
    # ---------------------------------------------------------

    @testset "divergence" begin
        expr = [x^2, y^2, z^2]
        
        div = divergence(expr, (x, y, z))
        divF =SymbolicFunction(div, (x, y, z))
        
        expected = 2x + 2y + 2z

        @test isequal(divF.expression, expected)
        @test divF(1.0, 2.0, 3.0) ≈ 12.0
    end

        # ---------------------------------------------------------
    # Mixed derivatives
    # ---------------------------------------------------------

    @testset "mixed derivatives" begin

        expr = x^2 * y + x * z^2 + y * z

        Hf = hessian(expr, (x, y, z))

        H = SymbolicFunction(Hf, (x, y, z))

        # d²f / dxdy = 2x
        @test isequal(simplify(H.expression[1,2] - 2x), 0)
        # d²f / dxdz = 2z
        @test isequal(simplify(H.expression[1,3] - 2z), 0)
        # d²f / dydz = 1
        @test isequal(simplify(H.expression[2,3] - 1), 0)

        result = H(2.0, 3.0, 4.0)

        expected = [
            2*3   2*2   2*4
            2*2   0     1
            2*4   1     2*2
        ]

        @test result ≈ Float64.(expected)
    end

    # ---------------------------------------------------------
    # Constant expression
    # ---------------------------------------------------------

    @testset "constant" begin

        expr = 42

        f = SymbolicFunction(expr, (x, y, z))

        @test f(1.0, 2.0, 3.0) ≈ 42.0
        @test f(-10.0, 100.0, 7.0) ≈ 42.0

        lf = laplacian(f.expression, f.variables)
        Δf = SymbolicFunction(lf, f.variables)

        @test Δf(1.0, 2.0, 3.0) ≈ 0.0
    end

    # ---------------------------------------------------------
    # More than three variables
    # ---------------------------------------------------------

    @testset "n-dimensional function" begin

        @variables x1 x2 x3 x4 x5

        variables = (x1, x2, x3, x4, x5)

        expr = x1^2 + x2^2 + x3^2 + x4^2 + x5^2

        f = SymbolicFunction(expr, variables)

        @test f(1, 2, 3, 4, 5) ≈ 55

        g = gradient(f.expression, f.variables)
        h = SymbolicFunction(g, f.variables)

        @test h(1, 2, 3, 4, 5) ≈ [2, 4, 6, 8, 10]

        lf = laplacian(f.expression, f.variables)
        Δf = SymbolicFunction(lf, f.variables)

        @test Δf(1, 2, 3, 4, 5) ≈ 10
    end
end


@testset "Symbolic Differential Operators (Symbolic Function Dispatch)" begin

    @variables x y z

    expr = x^2 + y^2 + z^2

    f = SymbolicFunction(expr, (x, y, z))


    # ========================================================
    # gradient
    # ========================================================

    @testset "gradient" begin

        g = gradient(f)

        @test g isa SymbolicFunction
        @test g.variables == f.variables

        @test isequal(
            g.expression,
            [2x, 2y, 2z]
        )

        @test g(1.0, 2.0, 3.0) ≈
              [2.0, 4.0, 6.0]
    end


    # ========================================================
    # laplacian
    # ========================================================

    @testset "laplacian" begin

        Δ = laplacian(f)

        @test Δ isa SymbolicFunction
        @test Δ.variables == f.variables

        @test isequal(
            Δ.expression,
            6
        )

        @test Δ(1.0, 2.0, 3.0) ≈ 6.0
        @test Δ(-10.0, 4.0, 7.0) ≈ 6.0
    end


    # ========================================================
    # hessian
    # ========================================================

    @testset "hessian" begin

        H = hessian(f)

        @test H isa SymbolicFunction
        @test H.variables == f.variables

        @test isequal(
            H.expression,
            [
                2 0 0
                0 2 0
                0 0 2
            ]
        )

        @test H(1.0, 2.0, 3.0) == [
            2.0 0.0 0.0
            0.0 2.0 0.0
            0.0 0.0 2.0
        ]
    end


    # ========================================================
    # Explicit variables
    # ========================================================

    @testset "explicit variables" begin

        variables = (x, y, z)

        g = gradient(f, variables)
        Δ = laplacian(f, variables)
        H = hessian(f, variables)

        @test g isa SymbolicFunction
        @test Δ isa SymbolicFunction
        @test H isa SymbolicFunction

        @test isequal(
            g.expression,
            [2x, 2y, 2z]
        )

        @test isequal(
            Δ.expression,
            6
        )

        @test isequal(
            H.expression,
            [
                2 0 0
                0 2 0
                0 0 2
            ]
        )
    end


    # ========================================================
    # Operator composition
    # ========================================================

    @testset "operator composition" begin

        g = gradient(f)
        H = hessian(f)
        Δ = laplacian(f)

        @test g isa SymbolicFunction
        @test H isa SymbolicFunction
        @test Δ isa SymbolicFunction

        @test isequal(g.expression, [2x, 2y, 2z])

        @test isequal(
            H.expression,
            [
                2 0 0
                0 2 0
                0 0 2
            ]
        )

        @test isequal(Δ.expression, 6)

        @test g(1.0, 2.0, 3.0) ≈ [2.0, 4.0, 6.0]

        @test H(1.0, 2.0, 3.0) ≈ [
            2.0 0.0 0.0
            0.0 2.0 0.0
            0.0 0.0 2.0
        ]

        @test Δ(1.0, 2.0, 3.0) ≈ 6.0
    end
    @testset "operator chaining" begin

        f = SymbolicFunction(
            x^3 + y^3 + z^3,
            (x, y, z)
        )

        g = gradient(f)
        H = hessian(f)
        Δ = laplacian(f)

        @test g isa SymbolicFunction
        @test H isa SymbolicFunction
        @test Δ isa SymbolicFunction

        # ∇(Δf) = ∇(6x + 6y + 6z)
        ∇Δ = gradient(Δ)

        @test ∇Δ isa SymbolicFunction

        @test isequal(
            ∇Δ.expression,
            [6, 6, 6]
        )

        @test ∇Δ(1.0, 2.0, 3.0) ≈ [6.0, 6.0, 6.0]
    end
end


