using SimplexGridFactory
using ExtendableGrids
using LinearAlgebra
using TetGen

function tetrahedralization_of_cube()
    
    builder=SimplexGridBuilder(Generator=TetGen)

    p1=point!(builder,0,0,0)
    p2=point!(builder,1,0,0)
    p3=point!(builder,1,1,0)
    p4=point!(builder,0,1,0)
    p5=point!(builder,0,0,1)
    p6=point!(builder,1,0,1)
    p7=point!(builder,1,1,1)
    p8=point!(builder,0,1,1)

    facetregion!(builder,1)
    facet!(builder,p1 ,p2 ,p3 ,p4)  
    facetregion!(builder,2)
    facet!(builder,p5 ,p6 ,p7 ,p8)  
    facetregion!(builder,3)
    facet!(builder,p1 ,p2 ,p6 ,p5)  
    facetregion!(builder,4)
    facet!(builder,p2 ,p3 ,p7 ,p6)  
    facetregion!(builder,5)
    facet!(builder,p3 ,p4 ,p8 ,p7)  
    facetregion!(builder,6)
    facet!(builder,p4 ,p1 ,p5 ,p8)

    meshy = simplexgrid(builder,maxvolume=0.25)

    return meshy
end

meshy = tetrahedralization_of_cube()
xVec = meshy.components[Coordinates][1,:]
yVec = meshy.components[Coordinates][2,:]
zVec = meshy.components[Coordinates][3,:]

function ufunc(xVec,yVec,zVec)
    ufunct = ones(length(xVec),length(yVec),length(zVec))
    for i = 1:length(xVec)
        for j = 1:length(yVec)
            for k = 1:length(zVec)
                ufunct[i] = xVec[i]^2 + yVec[j]^2 + zVec[k]^2
            end
        end
    end
    return ufunct
end

function Li_xBasis(x,xx,xi)
    lxi = 1
    for j = 1:length(xx)
        if xi != xx[j]
            lxi = lxi*(x-xx[j])/(xi-xx[j])
        end
    end
    return lxi
end

function Li_yBasis(y,yy,yi)
    lyi = 1
    for j = 1:length(yy)
        if yi != yy[j]
            lyi = lyi*(y-yy[j])/(yi-yy[j])
        end
    end
    return lyi
end

function Li_zBasis(z,zz,zi)
    lzi = 1
    for j = 1:length(zz)
        if zi != zz[j]
            lzi = lzi*(z-zz[j])/(zi-zz[j])
        end
    end
    return lzi
end

function usoln(x,y,xi,yi,xx,yy,zz,uij)
    lxi = Li_xBasis(x,xx,xi)
    lyi = Li_yBasis(y,yy,yi)
    lzi = Li_zBasis(z,zz,zi)
    u_soln = uij*lxi*lyi*lzi
    return u_soln
end

function usum(x,y,z,xx,yy,zz,uGrid)
    sum = 0
    for i = 1:length(xx)
        for j=1:length(yy)
            for k =1:length(zz)
                lzi = Li_zBasis(z,zz,zz[k])
                lyi = Li_yBasis(y,yy,yy[j])
                lxi = Li_xBasis(x,xx,xx[i])
                sum = sum + uGrid[i,j,k]*lxi*lyi*lzi
            end
        end
    end
    return sum
end

function unitPlot(xVec,yVec,zVec,xx,yy,zz,uGrid)
    plotty = zeros(length(xVec))
    for i = 1:length(xVec)
        plotty[i] = usum(xVec[i],yVec[i],zVec[i],xx,yy,zz,uGrid) 
    end
    return plotty
end

xx = LinRange(0,1,100)
yy = LinRange(0,1,100)
zz = LinRange(0,1,100)
uGrid = ufunc(xx,yy,zz)
plotty = unitPlot(xVec,yVec,zVec,xx,yy,zz,uGrid)

using GridVisualize
using GLMakie

scalarplot(meshy, plotty, Plotter=GLMakie,flevel=0.5)