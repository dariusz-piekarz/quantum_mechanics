using Test
using QuantumMechanics
using Symbolics


@testset "KineticOperator" begin

    @testset "Constructor" begin

        T̂32 = KineticOperator(Float32)
        T̂64 = KineticOperator(Float64)

        @test T̂32 isa KineticOperator{Float32}
        @test T̂64 isa KineticOperator{Float64}

        @test T̂32.coefficient == -0.5f0
        @test T̂64.coefficient == -0.5
    end


    @testset "1D function" begin

        @variables x

        ψ = SymbolicFunction(x^2, (x,))
        T̂ = KineticOperator(Float64)

        Tψ = T̂(ψ)

        @test Tψ isa SymbolicFunction
        @test Tψ.variables == ψ.variables

        # d²/dx² x² = 2
        # Tψ = -1/2 * 2 = -1
        @test Tψ.function_object(0.0) ≈ -1.0
        @test Tψ.function_object(1.0) ≈ -1.0
        @test Tψ.function_object(10.0) ≈ -1.0
    end


    @testset "3D function" begin

        @variables x y z

        # ψ = x² + y² + z²
        ψ = SymbolicFunction(
            x^2 + y^2 + z^2,
            (x, y, z)
        )

        T̂ = KineticOperator(Float64)

        Tψ = T̂(ψ)

        # ∇²ψ = 2 + 2 + 2 = 6
        # Tψ = -1/2 * 6 = -3

        @test Tψ.function_object(0.0, 0.0, 0.0) ≈ -3.0
        @test Tψ.function_object(1.0, 2.0, 3.0) ≈ -3.0
        @test Tψ.function_object(-10.0, 5.0, 2.0) ≈ -3.0
    end


    @testset "Non-trivial function" begin

        @variables x y z

        ψ = SymbolicFunction(
            x^2 * y + sin(z),
            (x, y, z)
        )

        T̂ = KineticOperator(Float64)

        Tψ = T̂(ψ)

        # ∇²(x²y + sin(z))
        #
        # d²/dx² (x²y) = 2y
        # d²/dy² (x²y) = 0
        # d²/dz² sin(z) = -sin(z)
        #
        # ∇²ψ = 2y - sin(z)
        #
        # Tψ = -1/2 (2y - sin(z))
        #     = -y + 1/2 sin(z)

        @test Tψ.function_object(0.0, 0.0, 0.0) ≈ 0.0

        @test Tψ.function_object(
            0.0,
            2.0,
            0.0
        ) ≈ -2.0

        @test Tψ.function_object(
            0.0,
            2.0,
            π / 2
        ) ≈ -1.5
    end


    @testset "Float32" begin

        @variables x y z

        ψ = SymbolicFunction(
            x^2 + y^2 + z^2,
            (x, y, z)
        )

        T̂ = KineticOperator(Float32)
        Tψ = T̂(ψ)

        @test T̂.coefficient isa Float32

        result = Tψ.function_object(
            1.0f0,
            2.0f0,
            3.0f0
        )

        @test result ≈ -3.0f0
    end


    @testset "Batch evaluation" begin

        @variables x y z

        ψ = SymbolicFunction(
            x^2 + y^2 + z^2,
            (x, y, z)
        )

        T̂ = KineticOperator(Float64)
        Tψ = T̂(ψ)

        X = [
            0.0 0.0 0.0
            1.0 0.0 0.0
            0.0 2.0 0.0
            0.0 0.0 3.0
        ]

        out = Vector{Float64}(undef, size(X, 1))

        evaluate!(Tψ, X, out)

        @test out ≈ fill(-3.0, 4)
    end

end