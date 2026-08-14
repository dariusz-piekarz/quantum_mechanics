using GLMakie
using Meshing, GeometryBasics
GLMakie.activate!()


function show_isosurface(
    f,
    h,
    ξ;
    color=(:crimson, 0.5),
    isoval=0
)

    algo = MarchingCubes(; iso=isoval)

    s = [h(x, y, z) for x in ξ, y in ξ, z in ξ]

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