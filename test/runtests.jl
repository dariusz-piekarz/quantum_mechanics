using Test

@testset "QuantumMechanics" begin
    for file in sort(readdir(@__DIR__))
        if file == "runtests.jl"
            continue
        end
        if endswith(file, ".jl")
            include(joinpath(@__DIR__, file))
        end
    end
end
