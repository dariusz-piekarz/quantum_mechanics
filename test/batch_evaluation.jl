using Test
using Symbolics
using QuantumMechanics


@testset "Batch evaluation" begin

    @variables x y z a b


    # ========================================================
    # Scalar functions
    # ========================================================

    @testset "scalar functions" begin

        # ----------------------------------------------------
        # 1 variable
        # ----------------------------------------------------

        @testset "1 variable" begin

            expr = x^2 + sin(x)

            f = SymbolicFunction(
                expr,
                (x,)
            )

            X = [
                1.0
                2.0
                3.0
            ] |> x -> reshape(x, :, 1)

            expected = [
                1.0^2 + sin(1.0),
                2.0^2 + sin(2.0),
                3.0^2 + sin(3.0)
            ]

            result = evaluate(f, X)

            @test result isa Vector
            @test size(result) == (3,)
            @test result ≈ expected


            out = Vector{Float64}(undef, 3)

            returned = evaluate!(
                f,
                X,
                out
            )

            @test returned === out
            @test out ≈ expected
        end


        # ----------------------------------------------------
        # 2 variables
        # ----------------------------------------------------

        @testset "2 variables" begin

            expr = x^2 + sin(y)

            f = SymbolicFunction(
                expr,
                (x, y)
            )

            X = [
                1.0 0.0
                2.0 1.0
                3.0 2.0
            ]

            expected = [
                1.0^2 + sin(0.0),
                2.0^2 + sin(1.0),
                3.0^2 + sin(2.0)
            ]

            result = evaluate(f, X)

            @test size(result) == (3,)
            @test result ≈ expected


            out = Vector{Float64}(undef, 3)

            returned = evaluate!(
                f,
                X,
                out
            )

            @test returned === out
            @test out ≈ expected
        end


        # ----------------------------------------------------
        # 3 variables
        # ----------------------------------------------------

        @testset "3 variables" begin

            expr = x^2 + sin(y) + exp(z)

            f = SymbolicFunction(
                expr,
                (x, y, z)
            )

            X = [
                1.0 0.0 0.0
                2.0 1.0 0.0
                3.0 2.0 1.0
            ]

            expected = [
                1.0^2 + sin(0.0) + exp(0.0),
                2.0^2 + sin(1.0) + exp(0.0),
                3.0^2 + sin(2.0) + exp(1.0)
            ]

            result = evaluate(f, X)

            @test size(result) == (3,)
            @test result ≈ expected


            out = Vector{Float64}(undef, 3)

            returned = evaluate!(
                f,
                X,
                out
            )

            @test returned === out
            @test out ≈ expected
        end


        # ----------------------------------------------------
        # 4 variables
        # ----------------------------------------------------

        @testset "4 variables" begin

            expr = x + 2y + 3z + 4a

            f = SymbolicFunction(
                expr,
                (x, y, z, a)
            )

            X = [
                1.0 2.0 3.0 4.0
                2.0 3.0 4.0 5.0
                3.0 4.0 5.0 6.0
            ]

            expected = [
                1.0 + 2*2.0 + 3*3.0 + 4*4.0,
                2.0 + 2*3.0 + 3*4.0 + 4*5.0,
                3.0 + 2*4.0 + 3*5.0 + 4*6.0
            ]

            result = evaluate(f, X)

            @test size(result) == (3,)
            @test result ≈ expected


            out = Vector{Float64}(undef, 3)

            returned = evaluate!(
                f,
                X,
                out
            )

            @test returned === out
            @test out ≈ expected
        end


        # ----------------------------------------------------
        # Generic fallback: 5 variables
        # ----------------------------------------------------

        @testset "generic scalar fallback" begin

            expr = x + 2y + 3z + 4a + 5b

            f = SymbolicFunction(
                expr,
                (x, y, z, a, b)
            )

            X = [
                1.0 2.0 3.0 4.0 5.0
                2.0 3.0 4.0 5.0 6.0
                3.0 4.0 5.0 6.0 7.0
            ]

            expected = [
                1 + 4 + 9 + 16 + 25,
                2 + 6 + 12 + 20 + 30,
                3 + 8 + 15 + 24 + 35
            ]

            result = evaluate(f, X)

            @test size(result) == (3,)
            @test result ≈ expected


            out = Vector{Float64}(undef, 3)

            returned = evaluate!(
                f,
                X,
                out
            )

            @test returned === out
            @test out ≈ expected
        end
    end


    # ========================================================
    # Gradient
    # ========================================================

    @testset "gradient" begin

        expr = x^2 + y^2 + z^2

        f = SymbolicFunction(
            expr,
            (x, y, z)
        )

        g = gradient(f)

        X = [
            1.0 2.0 3.0
            2.0 3.0 4.0
            3.0 4.0 5.0
        ]

        result = evaluate(g, X)

        expected = [
            2.0 4.0 6.0
            4.0 6.0 8.0
            6.0 8.0 10.0
        ]

        # Current layout:
        #
        # result[:, i] = gradient at point i
        #
        @test size(result) == (3, 3)

        for i in 1:3
            @test result[:, i] ≈ expected[i, :]
        end


        # evaluate!

        out = Matrix{Float64}(
            undef,
            3,
            3
        )

        returned = evaluate!(
            g,
            X,
            out
        )

        @test returned === out

        for i in 1:3
            @test out[:, i] ≈ expected[i, :]
        end
    end


    # ========================================================
    # Hessian
    # ========================================================

    @testset "hessian" begin

        expr = x^2 + y^2 + z^2

        f = SymbolicFunction(
            expr,
            (x, y, z)
        )

        H = hessian(f)

        X = [
            1.0 2.0 3.0
            2.0 3.0 4.0
            3.0 4.0 5.0
        ]

        expected = [
            2.0 0.0 0.0
            0.0 2.0 0.0
            0.0 0.0 2.0
        ]

        result = evaluate(H, X)

        @test result isa Array{Float64, 3}

        # Current layout:
        #
        # result[:, :, i]
        #
        # = Hessian at point i

        @test size(result) == (3, 3, 3)

        for i in 1:3
            @test result[:, :, i] ≈ expected
        end


        # evaluate!

        out = Array{Float64, 3}(
            undef,
            3,
            3,
            3
        )

        returned = evaluate!(
            H,
            X,
            out
        )

        @test returned === out

        for i in 1:3
            @test out[:, :, i] ≈ expected
        end
    end


    # ========================================================
    # Nontrivial Hessian
    # ========================================================

    @testset "nontrivial hessian" begin

        expr =
            x^2 * y +
            x * z^2 +
            y * z

        f = SymbolicFunction(
            expr,
            (x, y, z)
        )

        H = hessian(f)

        X = [
            1.0 2.0 3.0
            2.0 3.0 4.0
        ]

        result = evaluate(H, X)

        expected1 = [
            4.0 2.0 6.0
            2.0 0.0 1.0
            6.0 1.0 2.0
        ]

        expected2 = [
            6.0 4.0 8.0
            4.0 0.0 1.0
            8.0 1.0 4.0
        ]

        @test size(result) == (3, 3, 2)

        @test result[:, :, 1] ≈ expected1
        @test result[:, :, 2] ≈ expected2
    end


    # ========================================================
    # Single point
    # ========================================================

    @testset "single point" begin

        expr = x^2 + y^2 + z^2

        f = SymbolicFunction(
            expr,
            (x, y, z)
        )

        X = [
            1.0 2.0 3.0
        ]


        # ----------------------------------------------------
        # Scalar
        # ----------------------------------------------------

        result = evaluate(f, X)

        @test size(result) == (1,)
        @test result ≈ [14.0]


        # ----------------------------------------------------
        # Gradient
        # ----------------------------------------------------

        g = gradient(f)

        result_g = evaluate(g, X)

        @test size(result_g) == (3, 1)

        @test result_g[:, 1] ≈ [
            2.0,
            4.0,
            6.0
        ]


        # ----------------------------------------------------
        # Hessian
        # ----------------------------------------------------

        H = hessian(f)

        result_H = evaluate(H, X)

        @test size(result_H) == (3, 3, 1)

        @test result_H[:, :, 1] ≈ [
            2.0 0.0 0.0
            0.0 2.0 0.0
            0.0 0.0 2.0
        ]
    end


    # ========================================================
    # Many points
    # ========================================================

    @testset "many points" begin

        expr = x^2 + y^2 + z^2

        f = SymbolicFunction(
            expr,
            (x, y, z)
        )

        g = gradient(f)
        H = hessian(f)

        n = 1000

        X = rand(
            n,
            3
        )


        # ----------------------------------------------------
        # Scalar
        # ----------------------------------------------------

        result_f = evaluate(
            f,
            X
        )

        @test size(result_f) == (n,)

        for i in 1:n

            @test result_f[i] ≈
                sum(X[i, :] .^ 2)

        end


        # ----------------------------------------------------
        # Gradient
        # ----------------------------------------------------

        result_g = evaluate(
            g,
            X
        )

        @test size(result_g) == (3, n)

        for i in 1:n

            @test result_g[:, i] ≈
                2 .* X[i, :]

        end


        # ----------------------------------------------------
        # Hessian
        # ----------------------------------------------------

        result_H = evaluate(
            H,
            X
        )

        @test size(result_H) == (
            3,
            3,
            n
        )

        expected_H = [
            2.0 0.0 0.0
            0.0 2.0 0.0
            0.0 0.0 2.0
        ]

        for i in 1:n

            @test result_H[:, :, i] ≈
                expected_H

        end
    end


    # ========================================================
    # Preallocated output
    # ========================================================

    @testset "preallocated output" begin

        expr = x + y + z

        f = SymbolicFunction(
            expr,
            (x, y, z)
        )

        X = [
            1.0 2.0 3.0
            4.0 5.0 6.0
        ]

        out = Vector{Float64}(
            undef,
            2
        )

        returned = evaluate!(
            f,
            X,
            out
        )

        @test returned === out
        @test out ≈ [
            6.0,
            15.0
        ]
    end


    # ========================================================
    # Empty batch
    # ========================================================

    @testset "empty batch" begin

        expr = x + y + z

        f = SymbolicFunction(
            expr,
            (x, y, z)
        )

        X = Matrix{Float64}(
            undef,
            0,
            3
        )

        result = evaluate(
            f,
            X
        )

        @test size(result) == (0,)
    end


    # ========================================================
    # Output type
    # ========================================================

    @testset "output type" begin

        expr = x^2 + y^2

        f = SymbolicFunction(
            expr,
            (x, y)
        )

        X = Float32[
            1.0 2.0
            3.0 4.0
        ]

        result = evaluate(
            f,
            X,
            Float32
        )

        @test result isa Vector{Float32}

        @test result ≈ Float32[
            5.0,
            25.0
        ]
    end


    # ========================================================
    # 4-dimensional array output
    # ========================================================

    @testset "4 variables non-scalar output" begin

        expr = x^2 + y^2 + z^2 + a^2

        f = SymbolicFunction(
            expr,
            (x, y, z, a)
        )

        g = gradient(f)

        X = [
            1.0 2.0 3.0 4.0
            2.0 3.0 4.0 5.0
        ]

        result = evaluate(
            g,
            X
        )

        @test size(result) == (4, 2)

        @test result[:, 1] ≈ [
            2.0,
            4.0,
            6.0,
            8.0
        ]

        @test result[:, 2] ≈ [
            4.0,
            6.0,
            8.0,
            10.0
        ]
    end

end