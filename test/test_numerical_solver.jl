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

    L_forward = laplacian_1d(n, Val(:forward), h, Float64)
    L_backward = laplacian_1d(n, Val(:backward), h, Float64)
    L_central = laplacian_1d(n, Val(:central), h, Float64)

    @test size(L_forward) == (n, n)
    @test size(L_backward) == (n, n)
    @test size(L_central) == (n, n)

    @test issparse(L_forward)
    @test issparse(L_backward)
    @test issparse(L_central)

    @test eltype(L_forward) == Float64
    @test eltype(L_backward) == Float64
    @test eltype(L_central) == Float64

    @test Matrix(L_forward) == (1 / h^2) * [
        1   0   0   0   0
       -2   1   0   0   0
        1  -2   1   0   0
        0   1  -2   1   0
        0   0   1  -2   1
    ]

    @test Matrix(L_backward) == (1 / h^2) * [
        1  -2   1   0   0
        0   1  -2   1   0
        0   0   1  -2   1
        0   0   0   1  -2
        0   0   0   0   1
    ]

    @test Matrix(L_central) == (1 / h^2) * [
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

    L = laplacian(n, Val(:central), h, Float64, 3)

    @test size(L) == (n^3, n^3)
    @test issparse(L)
    @test eltype(L) == Float64
end


@testset "raw_laplacian_1d - invalid n" begin

    @test_throws DimensionMismatch raw_laplacian_1d(0, 0.1, Float64)
end


@testset "laplacian - invalid arguments" begin

    @test_throws DimensionMismatch laplacian(3, Val(:central), 0.1, Float64, 0)
    @test_throws MethodError laplacian(5, Val(:wrong), 0.1, Float64)
end