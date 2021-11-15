using SimplexGridFactory, ExtendableGrids
using LinearAlgebra
using TetGen
using BenchmarkTools
using Polyester
using GridVisualize, GLMakie
using GeometryBasics

function tetrahedralization_of_cube(; xlims=(-1,1), ylims=(-1,1), zlims=(-1,1), maxvolume=.1)
    
    builder=SimplexGridBuilder(Generator=TetGen)

    p1=point!(builder, xlims[1], ylims[1], zlims[1])
    p2=point!(builder, xlims[2], ylims[1], zlims[1])
    p3=point!(builder, xlims[2], ylims[2], zlims[1])
    p4=point!(builder, xlims[1], ylims[2], zlims[1])
    p5=point!(builder, xlims[1], ylims[1], zlims[2])
    p6=point!(builder, xlims[2], ylims[1], zlims[2])
    p7=point!(builder, xlims[2], ylims[2], zlims[2])
    p8=point!(builder, xlims[1], ylims[2], zlims[2])

    facetregion!(builder,1)
    facet!(builder, p1 ,p2, p3 ,p4)  
    facetregion!(builder,2)
    facet!(builder, p5 ,p6, p7 ,p8)  
    facetregion!(builder,3)
    facet!(builder, p1, p2, p6 ,p5)  
    facetregion!(builder,4)
    facet!(builder, p2, p3 ,p7 ,p6)  
    facetregion!(builder,5)
    facet!(builder, p3, p4 ,p8 ,p7)  
    facetregion!(builder,6)
    facet!(builder, p4, p1 ,p5 ,p8)

    meshy = simplexgrid(builder,maxvolume=maxvolume)

    return meshy
end

xlims=(-pi,pi)
ylims=(-pi,pi)
zlims=(-pi,pi)
meshy = tetrahedralization_of_cube(xlims=xlims, ylims=ylims, zlims=zlims, maxvolume=.1)
coord = meshy.components[Coordinates]
cellnodes = meshy.components[CellNodes]

xVec = meshy.components[Coordinates][1,:]
yVec = meshy.components[Coordinates][2,:]
zVec = meshy.components[Coordinates][3,:]
func = zeros(size(xVec))

for i = 1:length(xVec)
    func[i] = sqrt(xVec[i]^2+yVec[i]^2+zVec[i]^2)
end

#makeplanes(xyzmin,xyzmax,x,y,z)

xyzmin = [-pi;-pi;-pi]
xyzmax = [pi;pi;pi]
#planes = GridVisualize.makeplanes(xyzmin,xyzmax,xVec,yVec,zVec)
planes = []

flevels = [1;1.5]

tet_deets = GridVisualize.marching_tetrahedra(meshy,func,planes,flevels)
pts = tet_deets[1]
trngls = tet_deets[2]
fvals = tet_deets[3]

makie_triangles = Makie.to_triangles(hcat(getindex.(trngls,1), getindex.(trngls,2), getindex.(trngls,3)))
makie_pts = Makie.to_vertices(pts)
iso_mesh = GeometryBasics.normal_mesh(makie_pts, makie_triangles)

Makie.mesh(iso_mesh,color = :red, transparency = true)