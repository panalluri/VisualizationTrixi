using SimplexGridFactory, ExtendableGrids
using LinearAlgebra
using TetGen
using BenchmarkTools
using Polyester
using GridVisualize, GLMakie
using GeometryBasics
using Colors

function tetrahedralization_of_cube(; xlims=(-1,1), ylims=(-1,1), zlims=(-1,1), maxvolume=.01)
    
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

function ufunc(xVec,yVec,zVec)
    uGrid = zeros(length(xVec),length(yVec),length(zVec))
    for i = 1:length(xVec)
        for j = 1:length(yVec)
            for k = 1:length(zVec)
                # Taylor-Green Vortex 
                # uGridx = sin(xVec[i])*cos(yVec[j])*cos(zVec[k])
                # uGridy = -cos(xVec[i])*sin(yVec[j])*cos(zVec[k])
                # uGridz = 0
                # uGrid[i,j,k] = sqrt(uGridx^2+uGridy^2+uGridz^2)

                # Version w/ u=x, v=y, w=z
                uGridx = xVec[i]
                uGridy = yVec[j]
                uGridz = zVec[k]
                uGrid[i,j,k] = sqrt(uGridx^2+uGridy^2+uGridz^2)
            end
        end
    end
    return uGrid
end

function bary(xx,yy,zz)
    # TODO: avoid allocating these temporary arrays for additional speedup
    wx = zeros(1,length(xx))
    wy = zeros(1,length(yy))
    wz = zeros(1,length(zz))
    for j = 1:length(xx)
        wjx = 1
        for k = 1:length(xx)
            if xx[j] != xx[k]
                wjx = wjx * 1 / (xx[j]-xx[k])
            end
        end
        wx[j] = wjx
    end
    for j = 1:length(yy)
        wjy = 1
        for k = 1:length(yy)
            if yy[j] != yy[k]
                wjy = wjy * 1 / (yy[j]-yy[k])
            end
        end
        wy[j] = wjy
    end
    for j = 1:length(zz)
        wjz = 1
        for k = 1:length(zz)
            if zz[j] != zz[k]
                wjz = wjz * 1 / (zz[j] - zz[k])
            end
        end
        wz[j] = wjz
    end
    return wx, wy, wz
end

function usum(wx,wy,wz,x,y,z,xx,yy,zz,uGrid)
    sum = 0
    sumid = 0
    sumjd = 0
    sumkd = 0

    for k =1:length(zz)
        sumkd = sumkd + wz[k] / (z - zz[k])
    end
    for j =1:length(yy)
        sumjd = sumjd + wy[j] / (y - yy[j])
    end
    for i =1:length(xx)
        sumid = sumid + wx[i] / (x - xx[i])
    end

    for i = 1:length(xx)
        for j = 1:length(yy)
            for k =1:length(zz)
                sum = sum + uGrid[i,j,k] * (wz[k]/(z - zz[k]))/sumkd * (wy[j]/(y - yy[j]))/sumjd * (wx[i]/(x - xx[i]))/sumid
            end
        end
    end

    return sum
end

function unitPlot(xVec,yVec,zVec,xx,yy,zz,uGrid)
    wx, wy, wz = bary(xx, yy, zz)
    plotty = zeros(length(xVec))
    @batch for i = 1:length(xVec)
        plotty[i] = usum(wx, wy, wz, xVec[i], yVec[i], zVec[i], xx, yy, zz, uGrid) 
    end
    return plotty
end

fmat = [1 1.15 1.25 1.5; 1 1.15 1.25 1.5; 1 1.15 1.25 1.5]
planes = []

xyzlims = ([-3*pi,-pi], [-pi,pi], [pi,3*pi]) # xyzlims[i]
list_of_flevels = [fmat[i,:] for i in 1:size(fmat,1)]

list_of_meshes=[]
list_of_colors=[]

for j = 1:3
    meshy = tetrahedralization_of_cube(xlims=xyzlims[j], ylims=xyzlims[j], zlims=xyzlims[j], maxvolume=1.5)
    coord = meshy.components[Coordinates]
    cellnodes = meshy.components[CellNodes]
    xVec,yVec,zVec = (meshy.components[Coordinates][i,:] for i = 1:3)

    # chebyshev nodes
    n = length(xVec)-1
    #n = 7
    θ = LinRange(xyzlims[j][1], xyzlims[j][2], n+1)
    xn = sort(cos.(θ))

    xx = xn .* xyzlims[j][1] # LinRange(xlims..., 15)
    yy = xn .* xyzlims[j][1] # LinRange(ylims..., 15)
    zz = xn .* xyzlims[j][1] # LinRange(zlims..., 15)
    uGrid = ufunc(xx,yy,zz)
    func = unitPlot(xVec,yVec,zVec,xx,yy,zz,uGrid)

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