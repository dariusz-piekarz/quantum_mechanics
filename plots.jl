using GLMakie
using LaTeXStrings


include(joinpath(@__DIR__, "wave_function.jl"))
include(joinpath(@__DIR__, "isosurface.jl"))

"""
    plot_legendre(l; num_points=150)

Plot the Legendre polynomial \$P_l(x)\$ over the interval \$[-1, 1]\$.

# Arguments
- `l::Integer`: degree of the Legendre polynomial.
- `num_points::Integer=150`: number of sampled points used for plotting.

# Returns
A Makie figure containing the plotted polynomial.
"""
function plot_legendre(l::Integer; num_points::Integer=150)
    x = range(-1, 1, length=num_points)
    P = legendre.(l, x)

    fig = Figure(size=(900, 800))
    ax = Axis(fig[1, 1], xlabel="x", ylabel=L"P_{%$l}(x)", title=L"Legendre Polynomial $P_{%$l}(x)$")
    lines!(ax, x, P, color=:blue, linewidth=2)
    fig
end


"""
    plot_associated_legendre(l, m; num_points=150)

Plot the associated Legendre polynomial \$P_l^m(x)\$ over the interval \$[-1, 1]\$.

# Arguments
- `l::Integer`: degree.
- `m::Integer`: order.
- `num_points::Integer=150`: number of sampled points used for plotting.

# Returns
A Makie figure containing the plotted associated Legendre polynomial.
"""
function plot_associated_legendre(l::Integer, m::Integer; num_points::Integer=150)
    x = range(-1, 1, length=num_points)
    P = associated_legendre.(l, m, x)

    fig = Figure(size=(900, 800))
    ax = Axis(fig[1, 1], xlabel="x", ylabel=L"P_{%$l}^{%$m}(x)", title=L"Associated Legendre Polynomial $P_{%$l}^{%$m}(x)$")
    lines!(ax, x, P, color=:blue, linewidth=2)
    fig
end


function plot_spherical_harmonic(
    l::Integer,
    m::Integer;
    num_points::Integer=150,
    draw::Symbol=:real
)

    # --------------------------------------------------------
    # Angular coordinates
    # --------------------------------------------------------

    θ = range(0, π, length=num_points)
    φ = range(0, 2π, length=num_points)

    Y = spherical_harmonic(l, m, θ, φ)

    # --------------------------------------------------------
    # Quantity
    # --------------------------------------------------------

    if draw == :real

        values = Float64.(real.(Y))
        title = L"Re $Y_{%$l}^{%$m}$"

    elseif draw == :imag

        values = Float64.(imag.(Y))
        title = L"Im $Y_{%$l}^{%$m}$"

    elseif draw == :abs

        values = Float64.(abs.(Y))
        title = L"|Y_{%$l}^{%$m}|"

    elseif draw == :probability

        values = Float64.(abs2.(Y))
        title = L"|Y_{%$l}^{%$m}|^2"

    else

        throw(
            ArgumentError(
                "draw must be :real, :imag, :abs or :probability"
            )
        )
    end

    # --------------------------------------------------------
    # Angular grid
    # --------------------------------------------------------

    Θ = reshape(θ, :, 1)
    Φ = reshape(φ, 1, :)

    # --------------------------------------------------------
    # Radius
    #
    # For an orbital shape we use |values|.
    # The sign is encoded by the color.
    # --------------------------------------------------------

    r = abs.(values)

    # --------------------------------------------------------
    # Cartesian coordinates
    # --------------------------------------------------------

    x = r .* sin.(Θ) .* cos.(Φ)
    y = r .* sin.(Θ) .* sin.(Φ)
    z = r .* cos.(Θ)

    # --------------------------------------------------------
    # Figure
    # --------------------------------------------------------

    fig = Figure(size=(900, 800))

    ax = Axis3(
        fig[1, 1],
        aspect=:data,
        xlabel="x",
        ylabel="y",
        zlabel="z",
        title=title
    )

    # --------------------------------------------------------
    # Color
    # --------------------------------------------------------

    maxval = maximum(abs.(values))

    surface!(
        ax,
        x,
        y,
        z,
        color=values,
        colormap=[:blue, :white, :red],
        colorrange=(-maxval, maxval)
    )

    fig
end


function plot_laguerre(n::Integer; num_points::Integer=150)
    x = range(-1, 1, length=num_points)
    P = laguerre.(n, x)

    fig = Figure(size=(900, 800))
    ax = Axis(fig[1, 1], xlabel="x", ylabel=L"P_{%$n}(x)", title=L"Laguerre Polynomial $P_{%$n}(x)$")
    lines!(ax, x, P, color=:blue, linewidth=2)
    fig
end


function plot_generalized_laguerre(n::Integer, α::Union{Num,Number}; num_points::Integer=150)
    x = range(-1, 1, length=num_points)
    P = generalized_lagguerre.(n, α, x)

    fig = Figure(size=(900, 800))
    ax = Axis(fig[1, 1], xlabel="x", ylabel=L"P_{%$n}^{%$α}(x)", title=L"Generalized Laguerre Polynomial $P_{%$n}^{%$α}(x)$")
    lines!(ax, x, P, color=:blue, linewidth=2)
    fig
end


function plot_orbital(
    n::Integer,
    l::Integer,
    m::Integer;
    extent::Float64=15.0,
    step::Float64=0.15,
    level::Float64=0.1
)

    function density(x, y, z)

        r̃ = sqrt(x^2 + y^2 + z^2)

        if r̃ == 0
            θ = 0.0
            φ = 0.0
        else
            θ = acos(clamp(z / r̃, -1.0, 1.0))
            φ = atan(y, x)
        end

        r = α₀ * r̃

        ψ = wave_function(n, l, m, r, θ, φ)

        return abs2(ψ)
    end

    ξ = (-extent):step:extent

    ρmax = maximum(
        density(x, y, z)
        for x in ξ, y in ξ, z in ξ
    )

    isoval = level * ρmax

    fig = Figure(size=(900, 900))

    show_isosurface(
        fig[1, 1],
        density,
        ξ;
        isoval=isoval
    )

    return fig
end