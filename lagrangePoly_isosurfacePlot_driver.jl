#specify packages to be used for code
using SimplexGridFactory, ExtendableGrids
using LinearAlgebra
using TetGen
using BenchmarkTools
using Polyester
using GridVisualize, GLMakie
using GeometryBasics
using Colors

#create geometry for domain of interest
function tetrahedralization_of_cube(; xlims=(-1,1), ylims=(-1,1), zlims=(-1,1), maxvolume=.01)
    #set default parameters for a unit cube w/ max element volume size of 0.01
    #these values can be varied when calling the function
    
    #use TetGen to build geometry
    builder=SimplexGridBuilder(Generator=TetGen)

    #create pts
    p1=point!(builder, xlims[1], ylims[1], zlims[1])
    p2=point!(builder, xlims[2], ylims[1], zlims[1])
    p3=point!(builder, xlims[2], ylims[2], zlims[1])
    p4=point!(builder, xlims[1], ylims[2], zlims[1])
    p5=point!(builder, xlims[1], ylims[1], zlims[2])
    p6=point!(builder, xlims[2], ylims[1], zlims[2])
    p7=point!(builder, xlims[2], ylims[2], zlims[2])
    p8=point!(builder, xlims[1], ylims[2], zlims[2])

    #create faces
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

    #mesh over domain
    meshy = simplexgrid(builder,maxvolume=maxvolume)

    #return mesh for whole domain
    return meshy
end

#calculate function values at each of the chebyshev node points
function ufunc(xx,yy,zz)

    #pre-allocate matrix to hold chebyshev node point function evaluations
    uGrid = zeros(length(xx),length(yy),length(zz))

    #for our specific test function, we need to find the midpt of our conecentric spheres
    mid = mean(xx)

    #find function values at each of the chebyshev pts
    #plot concentric spheres a test function for the domain section of interest
    for i = 1:length(xx)
        for j = 1:length(yy)
            for k = 1:length(zz)
                uGridx = xx[i]-mid
                uGridy = yy[j]-mid
                uGridz = zz[k]-mid
                uGrid[i,j,k] = sqrt(uGridx^2+uGridy^2+uGridz^2)
            end
        end
    end

    #return matrix holding function evaluations of all chebyshev pts on a grid
    return uGrid
end

#find wj weight values in x-, y-, z-directions for barycentric lagrange polynomial implementation
#reference paper: https://people.maths.ox.ac.uk/trefethen/barycentric.pdf
function bary(xx,yy,zz)
    
    #pre-allocate vectors to hold weights in x-, y-, z-directions
    wx = zeros(1,length(xx))
    wy = zeros(1,length(yy))
    wz = zeros(1,length(zz))

    #calculate weights wj in x-direction
    for j = 1:length(xx)
        wjx = 1
        for k = 1:length(xx)
            #exclude the case where we divide by zero
            if !(xx[j] ≈ xx[k])
                wjx = wjx * 1 / (xx[j]-xx[k])
            end
        end
        wx[j] = wjx
    end

    #calculate weights wj in y-direction
    for j = 1:length(yy)
        wjy = 1
        for k = 1:length(yy)
            #exclude the case where we divide by zero
            if !(yy[j] ≈ yy[k])
                wjy = wjy * 1 / (yy[j]-yy[k])
            end
        end
        wy[j] = wjy
    end

    #calculate weights wj in z-direction
    for j = 1:length(zz)
        wjz = 1
        for k = 1:length(zz)
            #exclude the case where we divide by zero
            if !(zz[j] ≈ zz[k])
                wjz = wjz * 1 / (zz[j] - zz[k])
            end
        end
        wz[j] = wjz
    end

    #return vectors holding wj wight values in x-, y-, z-directions
    return wx, wy, wz
end

#calculate function values at the nodal pts found in the original TetGen 3D mesh generation
function unitPlot(xVec,yVec,zVec,xx,yy,zz,uGrid)

    #find the wj weight values in the x-, y-, z-directions
    wx, wy, wz = bary(xx, yy, zz)

    #pre-allocate a vector to hold function values at each of the TetGen nodal pts
    plotty = zeros(length(xVec))

    #find function values for TetGen nodal pts
    @batch for g = 1:length(xVec)

        #specify the x,y,z coordinates of the TetGen nodal pts for this loop
        x=xVec[g]
        y=yVec[g]
        z=zVec[g]

        #set various sum variables to zero for each loop
            #this ensures the summations in equation 4.2 from 
            #https://people.maths.ox.ac.uk/trefethen/barycentric.pdf
            #can be calculated
        sum = 0
        sumid = 0
        sumjd = 0
        sumkd = 0
    
        #sum over wj weights in the z-direction while dividing by z-grid pts
            #this becomes part of the denominator of eqn 4.2
        for k =1:length(zz)
            #avoid dividing by zero
            if !(z ≈ zz[k])
                sumkd = sumkd + wz[k] / (z - zz[k])
            end
        end

        #sum over wj weights in the y-direction while dividing by y-grid pts
            #this becomes part of the denominator of eqn 4.2
        for j =1:length(yy)
            #avoid dividing by zero
            if !(y ≈ yy[j])
                sumjd = sumjd + wy[j] / (y - yy[j])
            end
        end

        #sum over wj weights in the x-direction while dividing by x-grid pts
            #this becomes part of the denominator of eqn 4.2
        for i =1:length(xx)
            #avoid dividing by zero
            if !(x ≈ xx[i])
                sumid = sumid + wx[i] / (x - xx[i])
            end
        end
    
        #put numerator and denominator of eqn 4.2 together
            #multiply summations together as needed
        for i = 1:length(xx)
            for j = 1:length(yy)
                for k =1:length(zz)
                    #avoid dividing by zero
                    if !(z ≈ zz[k]) && !(y ≈ yy[j]) && !(x ≈ xx[i])
                        sum = sum + uGrid[i,j,k] * (wz[k]/(z - zz[k]))/sumkd * (wy[j]/(y - yy[j]))/sumjd * (wx[i]/(x - xx[i]))/sumid
                    end
                end
            end
        end

        #store the approximation for each TetGen nodal point in the plotty vector
        plotty[g] = sum
    end

    #return a vector containing function values for the TetGen nodal pts
    return plotty
end

#set limits for each of the unit cubes
xyzlims = ([-3*pi,-pi], [-pi,pi], [pi,3*pi]) 

#specify isosurfaces of interest for each unit block
fmat = [1 1.15 1.25 1.5; 1 1.15 1.25 1.5; 1 1.15 1.25 1.5]
list_of_flevels = [fmat[i,:] for i in 1:size(fmat,1)]

#no planar surfaces are of interest, other than the isosurfaces specified
    #this is important to note, as a non-empty planes vector will alter the 
    #output of the marching_tetrahedra algorithm
planes = []

#create empty vectors to hold isosurface mesh data and the colors corresponding to each mesh
list_of_meshes=[]
list_of_colors=[]

#find isosurfaces for each unit cube of interest
for j = 1:3

    #mesh over 3D domain
    meshy = tetrahedralization_of_cube(xlims=xyzlims[j], ylims=xyzlims[j], zlims=xyzlims[j], maxvolume=1)
    #extract coordinates and connectivity matrix from mesh
    coord = meshy.components[Coordinates]
    cellnodes = meshy.components[CellNodes] #connectivity matrix 
    #find point-wise coordinates on mesh split into x-, y-, z-components
    xVec,yVec,zVec = (meshy.components[Coordinates][i,:] for i = 1:3)

    #create chebyshev nodes
    n = length(xVec)-1
    n = 7
    θ = LinRange(0, pi, n+1)    
    xn = sort(cos.(θ))

    #set chebyshev nodes over interval of interest
    xdiff = 0.5 * (xyzlims[j][2] - xyzlims[j][1])
    xx = yy = zz = @. xn * xdiff + (xdiff + xyzlims[j][1])
    
    #calculate function values at chebyshev node points
    uGrid = ufunc(xx,yy,zz) 
    #find function values at node points from original TetGen mesh 
        #note that this is done using the chebyshev node grid and the function
        #values at these nodes
    func = unitPlot(xVec,yVec,zVec,xx,yy,zz,uGrid)

    #find each isosurface for a given unit block
    for i=1:length(list_of_flevels)

        #set isosurface value of interest
        level = list_of_flevels[i]
        
        #extract details about location isosurface pts using marching tetrahedra algorithm
        tet_deets = GridVisualize.marching_tetrahedra(meshy,func,planes,level)
        pts, trngls, fvals = tet_deets

        #output mesh that encompasses isosurface
        makie_triangles = Makie.to_triangles(hcat(getindex.(trngls,1), getindex.(trngls,2), getindex.(trngls,3)))
        makie_pts = Makie.to_vertices(pts)
        iso_mesh = GeometryBasics.normal_mesh(makie_pts, makie_triangles)

        #add newly found mesh to list of meshes
        push!(list_of_meshes, iso_mesh)
        #add a color to correspond w/ mesh of interest
        push!(list_of_colors, RGBA(0.1f0, 0.0f0, 0.0f0, 0.1f0))
    end
end

#plot isosurface meshes for each unit cube
fig, ax, plt = Makie.mesh(list_of_meshes[1],color = list_of_colors[1])
for i = 2:length(list_of_meshes)
    Makie.mesh!(list_of_meshes[i],color = list_of_colors[i])
end

#display final figure
fig