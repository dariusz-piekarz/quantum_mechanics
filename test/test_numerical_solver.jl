using Test
using SparseArrays
using KrylovKit

@testset "raw_laplacian" begin

    n = 5
    h = 0.1

    L = raw_laplacian(n, h, Float64)

    @test size(L) == (n, n + 2)
    @test eltype(L) == Float64
    @test issparse(L)

    # Współczynniki drugiej różnicy:
    # [1, -2, 1] / h^2
    @test L[1, 1] == 1 / h^2
    @test L[1, 2] == -2 / h^2
    @test L[1, 3] == 1 / h^2

    # Środkowy wiersz
    @test L[3, 3] == 1 / h^2
    @test L[3, 4] == -2 / h^2
    @test L[3, 5] == 1 / h^2

    # Ostatni wiersz
    @test L[n, n] == 1 / h^2
    @test L[n, n + 1] == -2 / h^2
    @test L[n, n + 2] == 1 / h^2

    # Poza trzema elementami w każdym wierszu powinno być zero
    for i in 1:n
        @test count(!iszero, L[i, :]) == 3
    end
end


@testset "laplacian - cut types" begin

    n = 5
    h = 0.1

    L_start   = laplacian(n, Val(:start),   h, Float64)
    L_end     = laplacian(n, Val(:end),     h, Float64)
    L_twotail = laplacian(n, Val(:twotail), h, Float64)

    # Wszystkie wersje powinny być n × n
    @test size(L_start)   == (n, n)
    @test size(L_end)     == (n, n)
    @test size(L_twotail) == (n, n)

    # Wszystkie powinny pozostać sparse
    @test issparse(L_start)
    @test issparse(L_end)
    @test issparse(L_twotail)

    # Typ
    @test eltype(L_start)   == Float64
    @test eltype(L_end)     == Float64
    @test eltype(L_twotail) == Float64

    # :start -> kolumny 3:(n+2)
    @test Matrix(L_start) == Matrix(
        (1 / h^2) .* [
            1  0  0  0  0
           -2  1  0  0  0
            1 -2  1  0  0
            0  1 -2  1  0
            0  0  1 -2  1
        ]
    )

    # :end -> kolumny 1:n
    @test Matrix(L_end) == Matrix(
        (1 / h^2) .* [
            1 -2  1  0  0
            0  1 -2  1  0
            0  0  1 -2  1
            0  0  0  1 -2
            0  0  0  0  1
        ]
    )

    # :twotail -> kolumny 2:(n+1)
    @test Matrix(L_twotail) == Matrix(
        (1 / h^2) .* [
           -2  1  0  0  0
            1 -2  1  0  0
            0  1 -2  1  0
            0  0  1 -2  1
            0  0  0  1 -2
        ]
    )
end


@testset "laplacian - different floating point types" begin

    n = 5
    h32 = Float32(0.1)
    h64 = Float64(0.1)

    L32 = raw_laplacian(n, h32, Float32)
    L64 = raw_laplacian(n, h64, Float64)

    @test eltype(L32) == Float32
    @test eltype(L64) == Float64

    @test size(L32) == (n, n + 2)
    @test size(L64) == (n, n + 2)
end


@testset "raw_laplacian - invalid n" begin

    @test_throws DimensionMismatch raw_laplacian(0, 0.1, Float64)
end


@testset "laplacian - invalid cut type" begin

    @test_throws MethodError laplacian(5, Val(:wrong), 0.1, Float64)
end


@testset "eigen_solver" begin

    n = 5
    h = 1.0

    L = laplacian(n, Val(:twotail), h, Float64)

    eigenvals, eigen_vects = eigen_solver(L, Float64)

    # Powinno zwrócić n wartości własnych
    @test length(eigenvals) == n

    # Powinno zwrócić n wektorów własnych
    @test length(eigen_vects) == n

    # Każdy wektor powinien mieć wymiar n
    @test all(length(v) == n for v in eigen_vects)

    # Typy
    @test eltype(eigenvals) == Float64
    @test all(eltype(v) == Float64 for v in eigen_vects)

    # Wszystkie wartości własne powinny być rzeczywiste
    @test all(isreal.(eigenvals))

    # Wartości własne powinny spełniać Av ≈ λv
    for (λ, v) in zip(eigenvals, eigen_vects)
        @test L * v ≈ λ * v atol=1e-10
    end
end