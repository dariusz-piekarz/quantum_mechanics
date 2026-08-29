using Test
using Symbolics
using QuantumMechanics


@testset "CoulombPotential" begin

    # ========================================================
    # Basic construction
    # ========================================================

    @testset "Construction" begin

        for T in (Float32, Float64)

            V = CoulombPotential{T}(
                -one(T),
                [one(T)],
                zeros(T, 1, 3)
            )

            @test V isa CoulombPotential{T}
            @test V.charge === -one(T)
            @test eltype(V.source_charges) === T
            @test eltype(V.source_positions) === T
        end
    end


    # ========================================================
    # Single source at origin
    #
    # V = q * Q / r
    # ========================================================

    @testset "Single source" begin

        for T in (Float32, Float64)

            V = CoulombPotential{T}(
                -one(T),
                [one(T)],
                zeros(T, 1, 3)
            )

            x = T(3)
            y = T(4)
            z = zero(T)

            # r = 5
            expected = -one(T) / T(5)

            @test V(x, y, z) ≈ expected
            @test typeof(V(x, y, z)) === T
        end
    end


    # ========================================================
    # Electron-electron
    #
    # (-1) * (-1) / r = +1/r
    # ========================================================

    @testset "Electron-electron potential" begin

        for T in (Float32, Float64)

            V = CoulombPotential{T}(
                -one(T),
                [-one(T)],
                zeros(T, 1, 3)
            )

            x = T(3)
            y = T(4)
            z = zero(T)

            expected = one(T) / T(5)

            @test V(x, y, z) ≈ expected
            @test typeof(V(x, y, z)) === T
        end
    end


    # ========================================================
    # Electron-nucleus
    #
    # (-1) * (+1) / r = -1/r
    # ========================================================

    @testset "Electron-nucleus potential" begin

        for T in (Float32, Float64)

            V = CoulombPotential{T}(
                -one(T),
                [one(T)],
                zeros(T, 1, 3)
            )

            expected = -one(T) / T(5)

            @test V(T(3), T(4), zero(T)) ≈ expected
        end
    end


    # ========================================================
    # Nucleus-electron
    #
    # (+1) * (-1) / r = -1/r
    # ========================================================

    @testset "Nucleus-electron potential" begin

        for T in (Float32, Float64)

            V = CoulombPotential{T}(
                one(T),
                [-one(T)],
                zeros(T, 1, 3)
            )

            expected = -one(T) / T(5)

            @test V(T(3), T(4), zero(T)) ≈ expected
        end
    end


    # ========================================================
    # Source translation
    # ========================================================

    @testset "Translated source" begin

        for T in (Float32, Float64)

            positions = T[
                1 2 3
            ]

            V = CoulombPotential{T}(
                -one(T),
                [one(T)],
                positions
            )

            # Point (4, 6, 3)
            # Distance from source (1,2,3):
            #
            # sqrt(3² + 4² + 0²) = 5

            expected = -one(T) / T(5)

            @test V(T(4), T(6), T(3)) ≈ expected
        end
    end


    # ========================================================
    # Multiple sources
    # ========================================================

    @testset "Multiple sources" begin

        for T in (Float32, Float64)

            positions = T[
                0 0 0
                3 0 0
            ]

            charges = T[
                1
                2
            ]

            V = CoulombPotential{T}(
                -one(T),
                charges,
                positions
            )

            # Point (0,4,0)
            #
            # r₁ = 4
            # r₂ = 5
            #
            # V = -1/4 - 2/5

            expected = -one(T) / T(4) - T(2) / T(5)

            @test V(zero(T), T(4), zero(T)) ≈ expected
        end
    end


    # ========================================================
    # Float32 type preservation
    # ========================================================

    @testset "Float32 type preservation" begin

        T = Float32

        V = CoulombPotential{T}(
            -one(T),
            [one(T)],
            zeros(T, 1, 3)
        )

        result = V(
            T(1),
            T(1),
            T(1)
        )

        @test result isa Float32
        @test typeof(result) === T
    end


    # ========================================================
    # Float64 type preservation
    # ========================================================

    @testset "Float64 type preservation" begin

        T = Float64

        V = CoulombPotential{T}(
            -one(T),
            [one(T)],
            zeros(T, 1, 3)
        )

        result = V(
            T(1),
            T(1),
            T(1)
        )

        @test result isa Float64
        @test typeof(result) === T
    end


    # ========================================================
    # Symbolic evaluation
    # ========================================================

    @testset "Symbolic evaluation" begin

        T = Float64

        V = CoulombPotential{T}(
            -one(T),
            [one(T)],
            zeros(T, 1, 3)
        )

        @variables x y z

        expr = V(x, y, z)

        @test expr isa Num

        # For a source at origin:
        #
        # V(x,y,z) = -1/sqrt(x²+y²+z²)

        expected = -one(T) / sqrt(x^2 + y^2 + z^2)

        @test Symbolics.simplify(expr - expected) == 0
    end


    # ========================================================
    # SymbolicFunction integration
    # ========================================================

    @testset "SymbolicFunction integration" begin

        T = Float32

        V = CoulombPotential{T}(
            -one(T),
            [one(T)],
            zeros(T, 1, 3)
        )

        @variables x y z

        f = SymbolicFunction(
            V(x, y, z),
            (x, y, z)
        )

        result = f(
            T(3),
            T(4),
            T(0)
        )

        @test result isa Float32
        @test result ≈ -T(1) / T(5)
    end
end


@testset "HarmonicPotential" begin

    @testset "Construction" begin
        k = 2.0
        V = HarmonicPotential(k)

        @test V isa HarmonicPotential{Float64}
        @test V.k == k
    end


    @testset "Numerical evaluation" begin
        k = 2.0
        V = HarmonicPotential(k)

        @test V(0.0, 0.0, 0.0) == 0.0

        @test V(1.0, 0.0, 0.0) == 1.0
        @test V(0.0, 1.0, 0.0) == 1.0
        @test V(0.0, 0.0, 1.0) == 1.0

        @test V(1.0, 2.0, 3.0) == 14.0
        @test V(-1.0, -2.0, -3.0) == 14.0
    end


    @testset "Float32 evaluation" begin
        k = Float32(2.0)
        V = HarmonicPotential(k)

        result = V(
            Float32(1.0),
            Float32(2.0),
            Float32(3.0)
        )

        @test result isa Float32
        @test result == Float32(14.0)
    end


    @testset "Symbolic evaluation" begin
        k = 2.0
        V = HarmonicPotential(k)

        @variables x y z

        result = V(x, y, z)

        expected = 1.0 * (x^2 + y^2 + z^2)

        @test isequal(result, expected)
    end


    @testset "Symmetry" begin
        V = HarmonicPotential(2.0)

        x, y, z = 1.2, -2.3, 3.4

        @test V(x, y, z) ≈ V(-x, y, z)
        @test V(x, y, z) ≈ V(x, -y, z)
        @test V(x, y, z) ≈ V(x, y, -z)
    end


    @testset "Scaling with k" begin
        V1 = HarmonicPotential(1.0)
        V2 = HarmonicPotential(2.0)

        x, y, z = 1.0, 2.0, 3.0

        @test V2(x, y, z) ≈ 2 * V1(x, y, z)
    end
end