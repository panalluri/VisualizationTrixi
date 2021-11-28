using SimplexGridFactory, ExtendableGrids
using LinearAlgebra
using TetGen
using BenchmarkTools
using Polyester
using GridVisualize, GLMakie
using GeometryBasics
using Colors

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

# xlims=(-3*pi,-pi,pi,3*pi)
# ylims=(-3*pi,-pi,pi,3*pi)
# zlims=(-3*pi,-pi,pi,3*pi)
fmat = [1 1.15 1.25 1.5; 1 1.15 1.25 1.5; 1 1.15 1.25 1.5]
planes = []

xyzlims = ([-3*pi,-pi], [-pi,pi], [pi,3*pi]) # xyzlims[i]
list_of_flevels = [fmat[i,:] for i in 1:size(fmat,1)]

list_of_meshes=[]
# list_of_colors=RGBA{Float32}[]
list_of_colors=[]

for j = 1:3
    meshy = tetrahedralization_of_cube(xlims=xyzlims[j], ylims=xyzlims[j], zlims=xyzlims[j], maxvolume=.1)
    coord = meshy.components[Coordinates]
    cellnodes = meshy.components[CellNodes]
    xVec,yVec,zVec = (meshy.components[Coordinates][i,:] for i = 1:3)
    func = zeros(size(xVec))
    for i = 1:length(xVec)
        func[i] = sqrt((xVec[i]-(2*j-4)*pi)^2+(yVec[i]-(2*j-4)*pi)^2^2+(zVec[i]-(2*j-4)*pi)^2^2)
    end

    for level in list_of_flevels[j]
        tet_deets = GridVisualize.marching_tetrahedra(meshy,func,planes,level)
        pts, trngls, fvals = tet_deets

        makie_triangles = Makie.to_triangles(hcat(getindex.(trngls,1), getindex.(trngls,2), getindex.(trngls,3)))
        makie_pts = Makie.to_vertices(pts)
        iso_mesh = GeometryBasics.normal_mesh(makie_pts, makie_triangles)

        push!(list_of_meshes, iso_mesh)
        push!(list_of_colors, RGBA(0.1f0, 0.0f0, 0.0f0, 0.1f0))
    end
end

# soln = Makie.mesh(iso_mesh,color = randn(length(pts)), transparency = true)
Makie.mesh(list_of_meshes[1],color = list_of_colors[1])

for i = 2:length(list_of_meshes)
    Makie.mesh!(list_of_meshes[i],color = list_of_colors[i])
end