using TetGen, GridVisualize
using LinearAlgebra
using Trixi
using StructArrays
using OrdinaryDiffEq

#create struct to hold node pts, func pts, connectivity matrix
struct PlotData3DTriangulated{DataType, NodeType, FaceNodeType, FaceDataType, VariableNames, PlottingTriangulation}
  x::NodeType # physical nodal coordinates, size (num_plotting_nodes x num_elements)
  y::NodeType
  z::NodeType
  data::DataType
  t::PlottingTriangulation
  x_face::FaceNodeType
  y_face::FaceNodeType
  z_face::FaceNodeType
  face_data::FaceDataType
  variable_names::VariableNames
end

#run trixi simulation and output trixi data
trixi_include("C:\\Users\\Prani\\.julia\\packages\\Trixi\\IAU6j\\examples\\dgmulti_3d\\elixir_euler_taylor_green_vortex.jl", tspan=(0 ,0.1), polydeg=7)

# get RefElemData, MeshData from Trixi
rd = solver.basis
md = mesh.md

rho, rho_v1, rho_v2, rho_v3, E = StructArrays.components(sol.u[end])
v1 = rho_v1./rho
v2 = rho_v2./rho
v3 = rho_v3./rho

@unpack x, y, z = md

using GLMakie, GeometryBasics, Colors
using StartUpDG

#create warped domain
# generate reference triangulation
input = TetGen.RawTetGenIO{Cdouble}(pointlist=vcat(transpose.(rd.rstp)...))
triangulation = tetrahedralize(input, "Q")
#connectivity matrix
connectivity = triangulation.tetrahedronlist

#find plotting pts from trixi mesh nodes
xp, yp, zp = (x -> rd.Vp * x).((x, y, z))

# r, s, t = reference coordinates. 
function derivative(u, coordinate, rd, md)
  @unpack Dr, Ds, Dt = rd
  @unpack rxJ, sxJ, txJ, ryJ, syJ, tyJ, rzJ, szJ, tzJ, J = md
  if coordinate==1
    # du/dx = du/dr * dr/dx + du/ds * ds/dx + du/dt * dt/dx
    return (rxJ .* (Dr * u) + sxJ .* (Ds * u) + txJ .* (Dt * u)) ./ J
  elseif coordinate==2
    return (ryJ .* (Dr * u) + syJ .* (Ds * u) + tyJ .* (Dt * u)) ./ J
  else #if coordinate==3
    return (rzJ .* (Dr * u) + szJ .* (Ds * u) + tzJ .* (Dt * u)) ./ J
  end
end

dudx = derivative(v1, 1, rd, md)
dudy = derivative(v1, 2, rd, md)
dudz = derivative(v1, 3, rd, md)
dvdx = derivative(v2, 1, rd, md)
dvdy = derivative(v2, 2, rd, md)
dvdz = derivative(v2, 3, rd, md)
dwdx = derivative(v3, 1, rd, md)
dwdy = derivative(v3, 2, rd, md)
dwdz = derivative(v3, 3, rd, md)

#find Q criteria
func_old = zeros(size(x))
omega = zeros(3,3)
S = zeros(3,3)
for j = 1:(size(x)[2])
    for i = 1:(size(x)[1])
        omega[1,2] = 0.5*(dudy[i,j]-dvdx[i,j])
        omega[1,3] = 0.5*(dudz[i,j]-dwdx[i,j])
        omega[2,1] = 0.5*(dvdx[i,j]-dudy[i,j])
        omega[2,3] = 0.5*(dvdz[i,j]-dwdx[i,j])
        omega[3,1] = 0.5*(dwdx[i,j]-dudz[i,j])
        omega[3,2] = 0.5*(dwdy[i,j]-dvdz[i,j])

        S[1,1] = dudx[i,j]
        S[1,2] = 0.5*(dudy[i,j]-dvdx[i,j])
        S[1,3] = 0.5*(dudz[i,j]-dwdx[i,j])
        S[2,1] = 0.5*(dvdx[i,j]-dudy[i,j])
        S[2,2] = dvdy[i,j]
        S[2,3] = 0.5*(dvdz[i,j]-dwdx[i,j])
        S[3,1] = 0.5*(dwdx[i,j]-dudz[i,j])
        S[3,2] = 0.5*(dwdy[i,j]-dvdz[i,j])
        S[3,3] = dwdz[i,j]

        Q = 0.5*(norm(omega)^2-norm(S)^2)
        func_old[i,j] = Q
    end
end

#reference element data
ref_elem = PlotData3DTriangulated(x,y,z,func_old,connectivity,[],[],[],[],[])

#interpolate func matrix
func = rd.Vp * func_old

#interpolated data
inter_elem = PlotData3DTriangulated(xp,yp,zp,func,connectivity,[],[],[],[],[])

#extract details about location isosurface pts using marching tetrahedra algorithm
level = [-.5] 

# create isosurface meshes to plot 
function global_plotting_triangulation_makie(inter_elem, level)

  xp = inter_elem.x
  yp = inter_elem.y
  zp = inter_elem.z
  func = inter_elem.data
  connectivity = inter_elem.t

  plotting_coordinates = zeros(3, size(xp, 1))
  num_elements = size(xp, 2)
  list_of_meshes = Vector{GeometryBasics.Mesh{3, Float32}}(undef, num_elements)
  sk = 1
  planes = []

  for e in 1:size(func, 2) # for each column
    plotting_coordinates[1, :] .= xp[:, e]
    plotting_coordinates[2, :] .= yp[:, e]
    plotting_coordinates[3, :] .= zp[:, e]
    pts, trngls, fvals = GridVisualize.marching_tetrahedra(plotting_coordinates,
                                                          connectivity,
                                                          func[:, e], planes, level)

    #output mesh that encompasses isosurface
    if length(pts) > 0
      makie_triangles = Makie.to_triangles(hcat((getindex.(trngls,i) for i in 1:3)...))
      iso_mesh = GeometryBasics.normal_mesh(Makie.to_vertices(pts), makie_triangles)

      # add newly found mesh to list of meshes
      list_of_meshes[sk] = iso_mesh
      sk += 1
    end    
  end

  plotting_mesh = merge(list_of_meshes[1:sk-1])
  return plotting_mesh

end

plotting_mesh = global_plotting_triangulation_makie(inter_elem, level)

solution_z = getindex.(plotting_mesh.position, 3)
Makie.mesh(plotting_mesh, color=solution_z, shading=false)