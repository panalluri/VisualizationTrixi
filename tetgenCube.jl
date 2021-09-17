using TetGen
using GeometryBasics
using GeometryBasics: Mesh, QuadFace

# Construct a cube out of Quads
points = Point{3, Float64}[
    (0.0, 0.0, 0.0), (1.0, 0.0, 0.0),
    (1.0, 1.0, 0.0), (0.0, 1.0, 0.0),
    (0.0, 0.0, 1.0), (1.0, 0.0, 1.0),
    (1.0, 1.0, 1.0), (0.0, 1.0, 1.0)
]

facets = QuadFace{Cint}[
    1:4,
    5:8,
    [1,5,6,2],
    [2,6,7,3],
    [3, 7, 8, 4],
    [4, 8, 5, 1]
]

mesh = Mesh(points, meta(facets))
result = tetrahedralize(mesh, "vpq1.414a0.1")

using GLMakie

GLMakie.mesh(normal_mesh(result), color=(:blue, 0.1), transparency=true)
GLMakie.wireframe!(result)

nodes = result.position