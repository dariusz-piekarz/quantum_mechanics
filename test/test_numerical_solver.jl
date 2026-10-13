using Test
using SparseArrays
using QuantumMechanics

@testset "raw_laplacian_1d" begin

    n = 5
    h = 0.1

    L = raw_laplacian_1d(n, h, Float64)

    @test size(L) == (n, n + 2)
    @test eltype(L) == Float64
    @test issparse(L)
    @test Matrix(L) == (1 / h^2) * [
        1  -2  1   0   0   0   0
        0   1  -2  1   0   0   0
        0   0   1  -2  1   0   0
        0   0   0   1  -2  1   0
        0   0   0   0   1  -2  1
    ]

    for i in 1:n
        @test count(!iszero, L[i, :]) == 3
    end
end


@testset "laplacian_1d - cut types" begin

    n = 5
    h = 0.1

    L = laplacian_1d(n, h, Float64)

    @test size(L) == (n, n)
    @test issparse(L)
    @test eltype(L) == Float64
    @test Matrix(L) == (1 / h^2) * [
       -2   1   0   0   0
        1  -2   1   0   0
        0   1  -2   1   0
        0   0   1  -2   1
        0   0   0   1  -2
    ]
end


@testset "laplacian - multidimensional" begin

    n = 3
    h = 0.1

    L = laplacian(n, h, Float64, 3)

    @test size(L) == (n^3, n^3)
    @test issparse(L)
    @test eltype(L) == Float64
end


@testset "raw_laplacian_1d - invalid n" begin
    @test_throws DimensionMismatch raw_laplacian_1d(0, 0.1, Float64)
    @test_throws ArgumentError raw_laplacian_1d(3, 0.0, Float64)
    @test_throws ArgumentError raw_laplacian_1d(3, Inf, Float64)
end
