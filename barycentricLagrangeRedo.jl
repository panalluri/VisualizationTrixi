using SimplexGridFactory, ExtendableGrids
using LinearAlgebra
using TetGen
using BenchmarkTools
using Polyester
using GridVisualize, GLMakie

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

xlims=(-pi,pi)
ylims=(-pi,pi)
zlims=(-pi,pi)
meshy = tetrahedralization_of_cube(xlims=xlims, ylims=ylims, zlims=zlims, maxvolume=.01)
xVec = meshy.components[Coordinates][1,:]
yVec = meshy.components[Coordinates][2,:]
zVec = meshy.components[Coordinates][3,:]

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

# TODO: compute curl from a given vector field instead of specifying it analytically
function uCurl(xVec,yVec,zVec)
    curly = zeros(length(xVec),length(yVec),length(zVec))
    for i = 1:length(xVec)
        for j = 1:length(yVec)
            for k = 1:length(zVec)
                # Taylor-Greene Vortex 
                # uGridx = sin(xVec[i])*cos(yVec[j])*cos(zVec[k])
                # uGridy = -cos(xVec[i])*sin(yVec[j])*cos(zVec[k])
                # uGridz = 0
                # uGrid[i,j,k] = sqrt(uGridx^2+uGridy^2+uGridz^2)

                # Curl version
                uGridx = -cos(xVec[i]) * sin(yVec[j]) * sin(zVec[k])
                uGridy = -sin(xVec[i]) * cos(yVec[j]) * sin(zVec[k])
                uGridz = sin(xVec[i]) * sin(yVec[j]) * cos(zVec[k]) + sin(xVec[i]) * sin(yVec[j]) * cos(zVec[k])
                curly[i,j,k] = sqrt(uGridx^2 + uGridy^2 + uGridz^2)

                # partial derivatives
                # dudx = cos(xVec[i])*cos(yVec[j])*cos(zVec[k])
                # dudy = -sin(xVec[i])*sin(yVec[j])*cos(zVec[k])
                # dudz = -sin(xVec[i])*cos(yVec[j])*sin(zVec[k])
                # dvdx = sin(xVec[i])*sin(yVec[j])*cos(zVec[k])
                # dvdy = -cos(xVec[i])*cos(yVec[j])*cos(zVec[k])
                # dvdz = cos(xVec[i])*sin(yVec[j])*sin(zVec[k])
                # dwdx = 0
                # dwdy = 0
                # dwdz = 0

                # omega = zeros(3,3)
                # omega[1,2] = 0.5*(dudy-dvdx)
                # omega[1,3] = 0.5*(dudz-dwdx)
                # omega[2,1] = 0.5*(dvdx-dudy)
                # omega[2,3] = 0.5*(dvdz-dwdx)
                # omega[3,1] = 0.5*(dwdx-dudz)
                # omega[3,2] = 0.5*(dwdy-dvdz)

                # S = zeros(3,3)
                # S[1,1] = dudx
                # S[1,2] = 0.5*(dudy-dvdx)
                # S[1,3] = 0.5*(dudz-dwdx)
                # S[2,1] = 0.5*(dvdx-dudy)
                # S[2,2] = dvdy
                # S[2,3] = 0.5*(dvdz-dwdx)
                # S[3,1] = 0.5*(dwdx-dudz)
                # S[3,2] = 0.5*(dwdy-dvdz)
                # S[3,3] = dwdz

                # Q = 0.5*(norm(omega)^2-norm(S)^2)
                # curly[i,j,k] = Q
                
            end
        end
    end
    return curly
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
                # lix = wx[i]/(x - xx[i])
                # ljy = wy[j]/(y - yy[j])
                # lkz = wz[k]/(z - zz[k])/sumkd
                sum = sum + uGrid[i,j,k] * wz[k]/(z - zz[k])/sumkd * wy[j]/(y - yy[j])/sumjd * wx[i]/(x - xx[i])/sumid
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

# chebyshev nodes
n = 7
θ = LinRange(0, pi, n+1)
xn = sort(cos.(θ))

xx = xn .* xlims[1] # LinRange(xlims..., 15)
yy = xn .* ylims[1] # LinRange(ylims..., 15)
zz = xn .* zlims[1] # LinRange(zlims..., 15)
uGrid = ufunc(xx,yy,zz)
curly = uCurl(xx,yy,zz)
plotty = unitPlot(xVec,yVec,zVec,xx,yy,zz,uGrid)
curlyPlot = unitPlot(xVec,yVec,zVec,xx,yy,zz,curly)

scalarplot(meshy, plotty, Plotter=GLMakie,outlinealpha=0,levels=[1;1.5;2])
# scalarplot(meshy, plotty, Plotter=GLMakie, outlinealpha=0, flevel=1)