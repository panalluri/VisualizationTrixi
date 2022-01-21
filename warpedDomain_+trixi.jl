using TetGen, GridVisualize 
using LinearAlgebra
using Trixi
using StructArrays
using OrdinaryDiffEq

#run trixi simulation and output trixi data
trixi_include("C:\\Users\\Prani\\.julia\\packages\\Trixi\\IAU6j\\examples\\dgmulti_3d\\elixir_euler_taylor_green_vortex.jl",tspan=(0.0,0.1))
rho, rho_v1, rho_v2, rho_v3, E = StructArrays.components(sol.u[end])
v1 = rho_v1./rho
v2 = rho_v2./rho
v3 = rho_v3./rho
@unpack x, y, z = mesh.md

using GLMakie, GeometryBasics, Colors
using StartUpDG

#create warped domain
rd = RefElemData(Hex(), N=3)
# generate reference triangulation
input = TetGen.RawTetGenIO{Cdouble}(pointlist=vcat(transpose.(rd.rstp)...))
triangulation = tetrahedralize(input, "Q")
#connectivity matrix
connectivity = triangulation.tetrahedronlist

#find plotting pts from trixi mesh nodes
xp, yp, zp = (x -> rd.Vp * x).((x, y, z))

#find Q criteria
func_old = zeros(size(x))
for j = 1:(size(x)[2])
    for i = 1:(size(x)[1])
        dudx = cos(x[i,j])*cos(y[i,j])*cos(z[i,j])
        dudy = -sin(x[i,j])*sin(y[i,j])*cos(z[i,j])
        dudz = -sin(x[i,j])*cos(y[i,j])*sin(z[i,j])
        dvdx = sin(x[i,j])*sin(y[i,j])*cos(z[i,j])
        dvdy = -cos(x[i,j])*cos(y[i,j])*cos(z[i,j])
        dvdz = cos(x[i,j])*sin(y[i,j])*sin(z[i,j])
        dwdx = 0
        dwdy = 0
        dwdz = 0

        omega = zeros(3,3)
        omega[1,2] = 0.5*(dudy-dvdx)
        omega[1,3] = 0.5*(dudz-dwdx)
        omega[2,1] = 0.5*(dvdx-dudy)
        omega[2,3] = 0.5*(dvdz-dwdx)
        omega[3,1] = 0.5*(dwdx-dudz)
        omega[3,2] = 0.5*(dwdy-dvdz)

        S = zeros(3,3)
        S[1,1] = dudx
        S[1,2] = 0.5*(dudy-dvdx)
        S[1,3] = 0.5*(dudz-dwdx)
        S[2,1] = 0.5*(dvdx-dudy)
        S[2,2] = dvdy
        S[2,3] = 0.5*(dvdz-dwdx)
        S[3,1] = 0.5*(dwdx-dudz)
        S[3,2] = 0.5*(dwdy-dvdz)
        S[3,3] = dwdz

        Q = 0.5*(norm(omega)^2-norm(S)^2)
        func_old[i,j] = Q
    end
end

#interpolate func matrix
func = rd.Vp * func_old

# to plot
isosurface_color = func.^2

#extract details about location isosurface pts using marching tetrahedra algorithm
level = [-1 -.5 -.25 .25 .5 1]

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