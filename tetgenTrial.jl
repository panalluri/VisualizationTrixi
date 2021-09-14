#using TetGen
using GeometryBasics

rect = Rect2(0.0, 0.0, 1.0, 1.0)
m = GeometryBasics.mesh(rect)

using GLMakie

GLMakie.mesh(normal_mesh(m), color=(:blue, 0.1), transparency=true)
GLMakie.wireframe!(m)