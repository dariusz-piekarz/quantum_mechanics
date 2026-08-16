using GLMakie
using Meshing, GeometryBasics
GLMakie.activate!()

"""
    show_isosurface(f, s; color=(:crimson, 0.5), isoval=0.0)

Render a 3D isosurface for a scalar field `s` on a figure `f`.
"""
function show_isosurface(
    f,
    s::AbstractArray{<:Real,3};
    color=(:crimson, 0.5),
    isoval=0.0
)

    algo = MarchingCubes(; iso=isoval)

    vts, fcs = Meshing.isosurface(s, algo)

    mc = GeometryBasics.Mesh(
        Point3f.(vts),
        GeometryBasics.TriangleFace.(fcs)
    )

    return mesh(
        f,
        normal_mesh(mc);
        color=color,
        diffuse=Vec3f(0.8),
        specular=Vec3f(1.1),
        shininess=30f0,
        backlight=5f0,
        transparency=true,
        axis=(; show_axis=false)
    )
end

"""
    show_isosurface!(f, s; color=(:crimson, 0.5), isoval=0.0)

Add a 3D isosurface to an existing plot (mutating version).
"""
function show_isosurface!(
    f,
    s::AbstractArray{<:Real,3};
    color=(:crimson, 0.5),
    isoval=0.0
)

    algo = MarchingCubes(; iso=isoval)

    vts, fcs = Meshing.isosurface(s, algo)

    mc = GeometryBasics.Mesh(
        Point3f.(vts),
        GeometryBasics.TriangleFace.(fcs)
    )

    return mesh!(
        f,
        normal_mesh(mc);
        color=color,
        diffuse=Vec3f(0.8),
        specular=Vec3f(1.1),
        shininess=30f0,
        backlight=5f0,
        transparency=true
    )
end


"""
    show_isosurface(f, density, ξ; color=(:crimson, 0.5), isoval=0.0)

Evaluate a scalar density function on a Cartesian grid and visualize its contour.
"""
function show_isosurface(
    f,
    density::Function,
    ξ::AbstractVector;
    color=(:crimson, 0.5),
    isoval=0.0
)
    grid = collect(ξ)
    values = [density(x, y, z) for x in grid, y in grid, z in grid]
    return show_isosurface(f, values; color=color, isoval=isoval)
end


"""
    show_isosurface(f, density, ξ; color=(:crimson, 0.5), isoval=0.0)

Fallback wrapper for iterable grids that are not already `AbstractVector`s.
"""
function show_isosurface(
    f,
    density::Function,
    ξ;
    color=(:crimson, 0.5),
    isoval=0.0
)
    return show_isosurface(f, density, collect(ξ); color=color, isoval=isoval)
end