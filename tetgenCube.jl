using TetGen
using GeometryBasics
using GeometryBasics: Mesh, QuadFace

# Construct a cube out of Quads
points = Point{3, Float64}[
    (0.0, 0.0, 0.0), (1.0, 0.0, 0.0),
    (1.0, 1.0, 0.0), (0.0, 1.0, 0.0),
    (0.0, 0.0, 1.0), (1.0, 0.0, 1.0),
    (1.0, 1.0, 1.0), (0.0, 1.0, 1.0)
]

facets = QuadFace{Cint}[
    1:4,
    5:8,
    [1,5,6,2],
    [2,6,7,3],
    [3, 7, 8, 4],
    [4, 8, 5, 1]
]

mesh = Mesh(points, meta(facets))
#change value after -a to control volume of elements
result = tetrahedralize(mesh, "q1.1a0.01")

using GLMakie

GLMakie.mesh(normal_mesh(result), color=(:blue, 0.1), transparency=true)
GLMakie.wireframe!(result)

nodes = result.position

function nodeOrganizer(nodes)
    xVec=zeros(length(nodes))
    yVec=zeros(length(nodes))
    zVec=zeros(length(nodes))
    for i=1:length(nodes)
        xVec[i]=nodes[i][1]
        yVec[i]=nodes[i][2]
        zVec[i]=nodes[i][3]
    end
    return xVec, yVec, zVec
end

function ufunc(xVec,yVec,zVec)
    ufunct = ones(length(xVec),length(yVec),length(zVec))
    for i = 1:length(xVec)
        for j = 1:length(yVec)
            for k = 1:length(zVec)
                ufunct[i] = zVec[k]*cos(xVec[i]) + zVec[k]*sin(yVec[j])
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

xVec, yVec, zVec = nodeOrganizer(nodes)
xx = LinRange(0,1,100)
yy = LinRange(0,1,100)
zz = LinRange(0,1,100)
uGrid = ufunc(xx,yy,zz)
plotty = unitPlot(xVec,yVec,zVec,xx,yy,zz,uGrid)
meshscatter(xVec,yVec,zVec,color=plotty)