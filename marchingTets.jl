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

xlims=(-3*pi,-pi,pi,3*pi)
ylims=(-3*pi,-pi,pi,3*pi)
zlims=(-3*pi,-pi,pi,3*pi)
fmat = [1 1.15 1.25 1.5; 1 1.15 1.25 1.5; 1 1.15 1.25 1.5]
planes = []

meshy = tetrahedralization_of_cube(xlims=xlims[1:2], ylims=ylims[1:2], zlims=zlims[1:2], maxvolume=.1)
coord = meshy.components[Coordinates]
cellnodes = meshy.components[CellNodes]

xVec = meshy.components[Coordinates][1,:]
yVec = meshy.components[Coordinates][2,:]
zVec = meshy.components[Coordinates][3,:]
func = zeros(size(xVec))

for i = 1:length(xVec)
    func[i] = sqrt((xVec[i]+2*pi)^2+(yVec[i]+2*pi)^2+(zVec[i]+2*pi)^2)
end

flevels = (fmat[1,:])'

tet_deets = GridVisualize.marching_tetrahedra(meshy,func,planes,flevels[1])
pts = tet_deets[1]
trngls = tet_deets[2]
fvals = tet_deets[3]

makie_triangles = Makie.to_triangles(hcat(getindex.(trngls,1), getindex.(trngls,2), getindex.(trngls,3)))
makie_pts = Makie.to_vertices(pts)
iso_mesh = GeometryBasics.normal_mesh(makie_pts, makie_triangles)

# soln = Makie.mesh(iso_mesh,color = randn(length(pts)), transparency = true)
Makie.mesh(iso_mesh,color = RGBA(0.0,0.0,1.0,0.1))
#to get isosurface 2, run Makie.mesh!(iso_mesh,color = RGBA(0.0,0.0,1.0,0.1)) w/ new flevel vector

for i = 2:(length(flevels)-1)
    tet_deets2 = GridVisualize.marching_tetrahedra(meshy,func,planes,flevels[i])
    pts2 = tet_deets2[1]
    trngls2 = tet_deets2[2]
    fvals2 = tet_deets2[3]

    makie_triangles2 = Makie.to_triangles(hcat(getindex.(trngls2,1), getindex.(trngls2,2), getindex.(trngls2,3)))
    makie_pts2 = Makie.to_vertices(pts2)
    iso_mesh2 = GeometryBasics.normal_mesh(makie_pts2, makie_triangles2)

    Makie.mesh!(iso_mesh2,color = RGBA(0.0,0.5,0.0,0.1))
end

tet_deets = GridVisualize.marching_tetrahedra(meshy,func,planes,flevels[length(flevels)])
pts = tet_deets[1]
trngls = tet_deets[2]
fvals = tet_deets[3]

makie_triangles = Makie.to_triangles(hcat(getindex.(trngls,1), getindex.(trngls,2), getindex.(trngls,3)))
makie_pts = Makie.to_vertices(pts)
iso_mesh = GeometryBasics.normal_mesh(makie_pts, makie_triangles)

Makie.mesh!(iso_mesh,color = RGBA(0.1,0.0,0.0,0.1))

for j = 2:(length(xlims)-1)
    global meshyy = tetrahedralization_of_cube(xlims=xlims[j:j+1], ylims=ylims[j:j+1], zlims=zlims[j:j+1], maxvolume=.1)
    global coordy = meshyy.components[Coordinates]
    global cellnodesy = meshyy.components[CellNodes]

    global xVecy = meshyy.components[Coordinates][1,:]
    global yVecy = meshyy.components[Coordinates][2,:]
    global zVecy = meshyy.components[Coordinates][3,:]
    global funcy = zeros(size(xVecy))

    for i = 1:length(xVecy)
        funcy[i] = sqrt((xVecy[i]-(2*j-4)*pi)^2+(yVecy[i]-(2*j-4)*pi)^2^2+(zVecy[i]-(2*j-4)*pi)^2^2)
    end

    global flevelsy = (fmat[j,:])'

    global tet_deetsy = GridVisualize.marching_tetrahedra(meshyy,funcy,planes,flevelsy[1])
    global ptsy = tet_deetsy[1]
    global trnglsy = tet_deetsy[2]
    global fvalsy = tet_deetsy[3]

    global makie_trianglesy = Makie.to_triangles(hcat(getindex.(trnglsy,1), getindex.(trnglsy,2), getindex.(trnglsy,3)))
    global makie_ptsy = Makie.to_vertices(ptsy)
    global iso_meshy = GeometryBasics.normal_mesh(makie_ptsy, makie_trianglesy)

    # soln = Makie.mesh(iso_mesh,color = randn(length(pts)), transparency = true)
    Makie.mesh!(iso_meshy,color = RGBA(0.0,0.0,1.0,0.1))
    #to get isosurface 2, run Makie.mesh!(iso_mesh,color = RGBA(0.0,0.0,1.0,0.1)) w/ new flevel vector

    for i = 2:(length(flevels)-1)
        tet_deetsz = GridVisualize.marching_tetrahedra(meshyy,funcy,planes,flevelsy[i])
        ptsz = tet_deetsz[1]
        trnglsz = tet_deetsz[2]
        fvalsz = tet_deetsz[3]

        makie_trianglesz = Makie.to_triangles(hcat(getindex.(trnglsz,1), getindex.(trnglsz,2), getindex.(trnglsz,3)))
        makie_ptsz = Makie.to_vertices(ptsz)
        iso_meshz = GeometryBasics.normal_mesh(makie_ptsz, makie_trianglesz)

        Makie.mesh!(iso_meshz,color = RGBA(0.0,0.5,0.0,0.1))
    end

    tet_deetsy = GridVisualize.marching_tetrahedra(meshyy,funcy,planes,flevelsy[length(flevelsy)])
    ptsy = tet_deetsy[1]
    trnglsy = tet_deetsy[2]
    fvalsy = tet_deetsy[3]

    makie_trianglesy = Makie.to_triangles(hcat(getindex.(trnglsy,1), getindex.(trnglsy,2), getindex.(trnglsy,3)))
    makie_ptsy = Makie.to_vertices(ptsy)
    iso_meshy = GeometryBasics.normal_mesh(makie_ptsy, makie_trianglesy)

    Makie.mesh!(iso_meshy,color = RGBA(0.1,0.0,0.0,0.1))
end