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
result = tetrahedralize(mesh, "vpq1.414a0.1")

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

function Li_xBasis(x,xVec,xi)
    lxi = 1
    for j = 1:length(xVec)
        if xi != xVec[j]
            lxi = lxi*(x-xVec[j])/(xi-xVec[j])
        end
    end
    return lxi
end

function Li_yBasis(y,yVec,yi)
    lyi = 1
    for j = 1:length(yVec)
        if yi != yVec[j]
            lyi = lyi*(y-yVec[j])/(yi-yVec[j])
        end
    end
    return lyi
end

function Li_zBasis(z,zVec,zi)
    lzi = 1
    for j = 1:length(zVec)
        if zi != zVec[j]
            lzi = lzi*(z-zVec[j])/(zi-zVec[j])
        end
    end
    return lzi
end

function usoln(x,y,xi,yi,xVec,yVec,zVec,uij)
    lxi = Li_xBasis(x,xVec,xi)
    lyi = Li_yBasis(y,yVec,yi)
    lzi = Li_zBasis(z,zVec,zi)
    u_soln = uij*lxi*lyi*lzi
    return u_soln
end

function usum(x,y,z,xVec,yVec,zVec,uGrid)
    sum = 0
    for i = 1:length(xVec)
        for j=1:length(yVec)
            for k =1:length(zVec)
                lzi = Li_zBasis(z,zVec,zVec[k])
                lyi = Li_yBasis(y,yVec,yVec[j])
                lxi = Li_xBasis(x,xVec,xVec[i])
                sum = sum + uGrid[i,j,k]*lxi*lyi*lzi
            end
        end
    end
    return sum
end

function unitPlot(xVec,yVec,zVec,uGrid)
    plotty = zeros(length(xVec))
    for i = 1:length(xVec)
        plotty[i] = usum(xVec[i],yVec[i],zVec[i],xVec,yVec,zVec,uGrid) 
    end
    return plotty
end

xVec, yVec, zVec = nodeOrganizer(nodes)
xx = LinRange(0,1,100)
yy = LinRange(0,1,100)
zz = LinRange(0,1,100)
uGrid = ufunc(xx,yy,zz)
plotty = unitPlot(xVec,yVec,zVec,uGrid)
meshscatter(xVec,yVec,zVec,color=plotty)