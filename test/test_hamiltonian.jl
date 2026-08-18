using Test
using Symbolics
using QuantumMechanics


@testset "Hamiltonian" begin

    @variables x y z


    # ============================================================
    # Basic hydrogen Hamiltonian
    # ============================================================

    @testset "Hydrogen - Float64" begin

        T = Float64

        V = CoulombPotential(
            -one(T),
            [one(T)],
            zeros(T, 1, 3)
        )

        K = KineticOperator(T)

        ψ = SymbolicFunction(
            x^2 + y^2 + z^2,
            (x, y, z)
        )

        Hψ = Hamiltonian(V, K, ψ)

        @test Hψ isa SymbolicFunction
        @test Hψ.variables == ψ.variables

        # ψ = x² + y² + z²
        #
        # ∇²ψ = 6
        #
        # Tψ = -3
        #
        # V = -1/r
        #
        # Hψ = -3 - r

        @test Hψ.function_object(
            1.0, 0.0, 0.0
        ) ≈ -4.0

        @test Hψ.function_object(
            0.0, 1.0, 0.0
        ) ≈ -4.0

        @test Hψ.function_object(
            0.0, 0.0, 2.0
        ) ≈ -5.0
    end


    # ============================================================
    # Float32
    # ============================================================

    @testset "Hydrogen - Float32" begin

        T = Float32

        V = CoulombPotential(
            -one(T),
            [one(T)],
            zeros(T, 1, 3)
        )

        K = KineticOperator(T)

        ψ = SymbolicFunction(
            x^2 + y^2 + z^2,
            (x, y, z)
        )

        Hψ = Hamiltonian(V, K, ψ)

        @test Hψ isa SymbolicFunction

        result = Hψ.function_object(
            T(1),
            T(0),
            T(0)
        )

        @test result ≈ T(-4)
    end


    # ============================================================
    # Constant wavefunction
    # ============================================================

    @testset "Constant wavefunction" begin

        T = Float64

        V = CoulombPotential(
            -one(T),
            [one(T)],
            zeros(T, 1, 3)
        )

        K = KineticOperator(T)

        ψ = SymbolicFunction(
            1.0,
            (x, y, z)
        )

        Hψ = Hamiltonian(V, K, ψ)

        # ∇² 1 = 0
        #
        # Hψ = Vψ = -1/r

        @test Hψ.function_object(
            1.0, 0.0, 0.0
        ) ≈ -1.0

        @test Hψ.function_object(
            0.0, 2.0, 0.0
        ) ≈ -0.5
    end


    # ============================================================
    # Pure kinetic contribution
    # ============================================================

    @testset "Zero potential" begin

        T = Float64

        V = CoulombPotential(
            zero(T),
            [one(T)],
            zeros(T, 1, 3)
        )

        K = KineticOperator(T)

        ψ = SymbolicFunction(
            x^2 + y^2 + z^2,
            (x, y, z)
        )

        Hψ = Hamiltonian(V, K, ψ)

        # Hψ = -1/2 ∇²ψ
        #     = -3

        @test Hψ.function_object(
            0.0, 0.0, 0.0
        ) ≈ -3.0

        @test Hψ.function_object(
            10.0, 20.0, 30.0
        ) ≈ -3.0
    end


    # ============================================================
    # Multiple nuclei
    # ============================================================

    @testset "Multiple nuclei" begin

        T = Float64

        V = CoulombPotential(
            -one(T),
            [one(T), one(T)],
            [
                0.0 0.0 0.0
                2.0 0.0 0.0
            ]
        )

        K = KineticOperator(T)

        ψ = SymbolicFunction(
            1.0,
            (x, y, z)
        )

        Hψ = Hamiltonian(V, K, ψ)

        # At (1, 0, 0):
        #
        # r₁ = 1
        # r₂ = 1
        #
        # V = -2
        #
        # Tψ = 0
        #
        # Hψ = -2

        @test Hψ.function_object(
            1.0, 0.0, 0.0
        ) ≈ -2.0
    end


    # ============================================================
    # Batch evaluation
    # ============================================================

    @testset "Batch evaluation" begin

        T = Float64

        V = CoulombPotential(
            -one(T),
            [one(T)],
            zeros(T, 1, 3)
        )

        K = KineticOperator(T)

        ψ = SymbolicFunction(
            x^2 + y^2 + z^2,
            (x, y, z)
        )

        Hψ = Hamiltonian(V, K, ψ)

        X = [
            1.0 0.0 0.0
            0.0 1.0 0.0
            0.0 0.0 2.0
            3.0 0.0 0.0
        ]

        out = Vector{Float64}(undef, size(X, 1))

        evaluate!(Hψ, X, out)

        expected = [
            -4.0,
            -4.0,
            -5.0,
            -1.0 / 3.0 * 9.0 - 3.0
        ]

        @test out ≈ expected
    end

end