using GLMakie, GeometryBasics, Colors
using StartUpDG
using TetGen, GridVisualize 
using LinearAlgebra

#create warped domain
rd = RefElemData(Hex(), N=3)

# generate reference triangulation
input = TetGen.RawTetGenIO{Cdouble}(pointlist=vcat(transpose.(rd.rstp)...))
triangulation = tetrahedralize(input, "Q")

md = MeshData(uniform_mesh(Hex(), 2)..., rd) 

@unpack x, y, z = md 
xp, yp, zp = (x -> rd.Vp * x).((x, y, z))

# to plot
f(x, y, z) = sin(pi * x) * sin(pi * y) * sin(pi * z)
func = f.(xp, yp, zp)
isosurface_color = func.^2

#extract details about location isosurface pts using marching tetrahedra algorithm
level = [-.5 .5]

#create empty vectors to hold isosurface mesh data and the colors corresponding to each mesh
plotting_coordinates = zeros(3, size(xp, 1))
num_elements = size(xp, 2)
list_of_meshes = Vector{GeometryBasics.Mesh{3, Float32}}(undef, num_elements)
for e in 1:size(func, 2) # for each column 
    plotting_coordinates[1, :] .= xp[:, e]
    plotting_coordinates[2, :] .= yp[:, e]
    plotting_coordinates[3, :] .= zp[:, e]
    pts, trngls, fvals = GridVisualize.marching_tetrahedra(plotting_coordinates, 
                                                           triangulation.tetrahedronlist, 
                                                           func[:, e], [], level)

    #output mesh that encompasses isosurface
    makie_triangles = Makie.to_triangles(hcat((getindex.(trngls,i) for i in 1:3)...))    
    iso_mesh = GeometryBasics.normal_mesh(Makie.to_vertices(pts), makie_triangles)

    # add newly found mesh to list of meshes
    list_of_meshes[e] = iso_mesh
end
plotting_mesh = merge([list_of_meshes...])
solution_z = getindex.(plotting_mesh.position, 3)
Makie.mesh(plotting_mesh, color=solution_z)