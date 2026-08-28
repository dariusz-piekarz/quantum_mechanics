using Test
using Symbolics
using QuantumMechanics

# ============================================================
# Helpers
# ============================================================

const ATOL = 1e-10
const RTOL = 1e-8

# ============================================================
# Factorials
# ============================================================

@testset "Factorials" begin

    @test factorials_table(0) == [1]
    @test factorials_table(1) == [1, 1]
    @test factorials_table(5) == [1, 1, 2, 6, 24, 120]

    @test_throws ArgumentError factorials_table(-1)

    facts = factorials_table(10)

    @test resolve_factorials(facts, 5) === facts
    @test_throws ArgumentError resolve_factorials(facts, 11)
end


# ============================================================
# Legendre polynomials
# ============================================================

@testset "Legendre polynomials" begin

    # Basic values
    @test legendre(0, 0.5) ≈ 1.0
    @test legendre(1, 0.5) ≈ 0.5

    @test legendre(2, 0.5) ≈
          (3 * 0.5^2 - 1) / 2

    @test legendre(2, 0.0) ≈ -0.5

    # Known polynomials
    x = 0.37

    @test legendre(2, x) ≈
          (3x^2 - 1) / 2

    @test legendre(3, x) ≈
          (5x^3 - 3x) / 2

    @test legendre(4, x) ≈
          (35x^4 - 30x^2 + 3) / 8

    # Endpoints
    for l in 0:10
        @test legendre(l, 1.0) ≈ 1.0
        @test legendre(l, -1.0) ≈ (-1)^l
    end

    # Negative degree
    @test_throws ArgumentError legendre(-1, 0.5)

    # Array version
    xgrid = collect(range(-1.0, 1.0, length=101))

    P = legendre(5, xgrid)

    @test size(P) == size(xgrid)
    @test eltype(P) <: Number
    @test all(isfinite, P)

    # Scalar vs array implementation
    for i in eachindex(xgrid)
        @test P[i] ≈ legendre(5, xgrid[i])
    end
end


# ============================================================
# Associated Legendre polynomials
# ============================================================

@testset "Associated Legendre polynomials" begin

    # P_l^0(x) = P_l(x)
    x = 0.37

    for l in 0:8
        @test associated_legendre(l, 0, x) ≈
              legendre(l, x)
    end

    # Known values
    #
    # P_1^1(x) = -sqrt(1-x^2)
    #
    @test associated_legendre(1, 1, x) ≈
          -sqrt(1 - x^2)

    # P_2^1(x) = -3x sqrt(1-x^2)
    @test associated_legendre(2, 1, x) ≈
          -3x * sqrt(1 - x^2)

    # P_2^2(x) = 3(1-x^2)
    @test associated_legendre(2, 2, x) ≈
          3(1 - x^2)

    # Negative m relation
    #
    # P_l^{-m}(x)
    # =
    # (-1)^m (l-m)!/(l+m)! P_l^m(x)
    #
    facts = factorials_table(12)

    for l in 1:6
        for m in 1:l

            lhs = associated_legendre(l, -m, x, facts)

            rhs =
                (-1)^m *
                facts[l-m+1] /
                facts[l+m+1] *
                associated_legendre(l, m, x, facts)

            @test lhs ≈ rhs
        end
    end

    # Array implementation
    xgrid = collect(range(-1.0, 1.0, length=101))

    for l in 0:5
        for m in (-l):l

            P = associated_legendre(l, m, xgrid)

            @test size(P) == size(xgrid)
            @test all(isfinite, P)

            for i in eachindex(xgrid)
                @test P[i] ≈ associated_legendre(l, m, xgrid[i])
            end
        end
    end

    # Invalid quantum numbers
    @test_throws ArgumentError associated_legendre(-1, 0, 0.5)
    @test_throws ArgumentError associated_legendre(2, 3, 0.5)
end


# ============================================================
# Laguerre polynomials
# ============================================================

@testset "Laguerre polynomials" begin

    @test laguerre(0, 0.5) == 1.0
    @test laguerre(1, 0.5) ≈ 1 - 0.5

    @test laguerre(2, 0.5) ≈
          1 - 2 * 0.5 + 0.5^2 / 2

    @test laguerre(2, 2.0) ≈
          1 - 4.0 + 2.0^2 / 2

    x = 0.37

    # Explicit formulas
    @test laguerre(2, x) ≈
          1 - 2x + x^2 / 2

    @test laguerre(3, x) ≈
          1 - 3x + 3x^2 / 2 - x^3 / 6

    # Invalid degree
    @test_throws ArgumentError laguerre(-1, x)

    # Array implementation
    xgrid = collect(range(0.0, 10.0, length=101))

    L = laguerre.(5, xgrid)

    @test size(L) == size(xgrid)
    @test all(isfinite, L)

    # Scalar vs array
    for i in eachindex(xgrid)
        @test L[i] ≈ laguerre(5, xgrid[i])
    end
end


# ============================================================
# Generalized Laguerre polynomials
# ============================================================

@testset "Generalized Laguerre polynomials" begin

    @test generalized_laguerre(0, 0.5, 0.5) == 1.0

    @test generalized_laguerre(1, 0.5, 2.0) ≈
          1.5 - 2.0

    x = 0.7
    α = 2.0

    # L_0^(α)(x) = 1
    @test generalized_laguerre(0, α, x) ≈ 1.0

    # L_1^(α)(x) = α + 1 - x
    @test generalized_laguerre(1, α, x) ≈
          α + 1 - x

    # α = 0 reduces to ordinary Laguerre
    for n in 0:8
        @test generalized_laguerre(n, 0.0, x) ≈
              laguerre(n, x)
    end

    # Array implementation
    xgrid = collect(range(0.0, 10.0, length=101))

    L = generalized_laguerre(5, 2.0, xgrid)

    @test size(L) == size(xgrid)
    @test all(isfinite, L)

    # Invalid arguments
    @test_throws ArgumentError generalized_laguerre(-1, 1.0, x)
    @test_throws ArgumentError generalized_laguerre(2, -2.0, x)
end


# ============================================================
# Spherical harmonics
# ============================================================

@testset "Spherical harmonics" begin

    # Y_0^0 = 1 / sqrt(4π)
    @test spherical_harmonic(0, 0, 0.0, 0.0) ≈
          1 / sqrt(4π)

    # Normalization constant
    @test N(1, 0) ≈ sqrt(3 / (4π))

    θ = 0.7
    φ = 1.2

    # Y_l^0 should be real
    for l in 0:5
        Y = spherical_harmonic(l, 0, θ, φ)

        @test isreal(Y)
    end

    # |Y_l^m| is independent of φ
    for l in 0:5
        for m in (-l):l

            Y1 = spherical_harmonic(l, m, θ, 0.0)
            Y2 = spherical_harmonic(l, m, θ, 1.0)

            @test abs(Y1) ≈ abs(Y2)
        end
    end

    # Negative m relation:
    #
    # Y_l^{-m} = (-1)^m conj(Y_l^m)
    #
    for l in 1:5
        for m in 1:l

            lhs = spherical_harmonic(l, -m, θ, φ)

            rhs =
                (-1)^m *
                conj(spherical_harmonic(l, m, θ, φ))

            @test lhs ≈ rhs
        end
    end

    # Array implementation
    θgrid = collect(range(0.0, π, length=25))
    φgrid = collect(range(0.0, 2π, length=30))

    Y = spherical_harmonic(2, 1, θgrid, φgrid)

    @test size(Y) == (length(θgrid), length(φgrid))
    @test all(isfinite, Y)

    # Invalid quantum numbers
    @test_throws ArgumentError spherical_harmonic(-1, 0, θ, φ)
    @test_throws ArgumentError spherical_harmonic(2, 3, θ, φ)
end
