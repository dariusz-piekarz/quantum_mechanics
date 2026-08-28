using QuantumMechanics
using Symbolics

abstract type AbstractPotential end


struct CoulombPotential{T} <: AbstractPotential
    charge::T
    source_charges::Vector{T}
    source_positions::Matrix{T}
end


function (V::CoulombPotential)(x::T, y::T, z::T)::T where {T}

    result = zero(T)

    @inbounds for A in eachindex(V.source_charges)

        dx = x - V.source_positions[A, 1]
        dy = y - V.source_positions[A, 2]
        dz = z - V.source_positions[A, 3]

        r = sqrt(dx^2 + dy^2 + dz^2)

        result += V.charge * V.source_charges[A] / r
    end

    return result
end


function (V::CoulombPotential{T})(x::Num, y::Num, z::Num)::Num where {T}

    result = Num(0)

    @inbounds for A in eachindex(V.source_charges)

        dx = x - V.source_positions[A, 1]
        dy = y - V.source_positions[A, 2]
        dz = z - V.source_positions[A, 3]

        r = sqrt(dx^2 + dy^2 + dz^2)

        result += V.charge * V.source_charges[A] / r
    end

    return result
end


function V_e(::Type{T}) where {T}
    return CoulombPotential{T}(
        -one(T),
        [one(T)],
        zeros(T, 1, 3)
    )
end

function V_ee(::Type{T}) where {T}
    return CoulombPotential{T}(
        -one(T),
        [-one(T)],
        zeros(T, 1, 3)
    )
end

function V_en(::Type{T}) where {T}
    return CoulombPotential{T}(
        -one(T),
        [one(T)],
        zeros(T, 1, 3)
    )
end

function V_ne(::Type{T}) where {T}
    return CoulombPotential{T}(
        one(T),
        [-one(T)],
        zeros(T, 1, 3)
    )
end
