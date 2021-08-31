using Trixi
using StructArrays
using Plots

include("C:\\Users\\Prani\\.julia\\dev\\Trixi\\examples\\dgmulti_2d\\elixir_euler_kelvin_helmholtz_instability.jl")
rho, rho_v1, rho_v2, E = StructArrays.components(sol.u[end])
v1 = rho_v1./rho
v2 = rho_v2./rho
@unpack x, y = mesh.md;
x = mesh.md.x
y = mesh.md.y

function makeGrid(x,ele)
    xx = x[:,ele]
    xGrid = zeros(4,4)
    z = [4 3 2 1]
    for i = 1:4
        xGrid[z[i],:] = transpose(xx[(i-1)*4+1:i*4])
    end
    return xGrid
end

function v1_fit(x,y,v1,ele)
    sumFunc = 0
    xx = makeGrid(x,ele)
    yy = makeGrid(y,ele)
    v11 = makeGrid(v1,ele)
    xLine = xx[1,:]
    yLine = yy[:,1]
    #how do I store each polynomial below as an element in a vector?
    xPoly = zeros(4)
    yPoly = zeros(4)
    #Find x polynomials
    for i = 1:4
        sum1 = 0
        for j = 1:4
            #how do I use symbolic x?
            sum1 = sum1 + (x-xLine[j])/(xLine[i]-xLine[j])
        end
        xPoly[i] = sum1
    end
    #Find y polynomials
    for i = 1:4
        sum1 = 0
        for j = 1:4
            #how do I use symbolic y?
            sum1 = sum1 + (y-yLine[j])/(yLine[i]-yLine[j])
        end
        yPoly[i] = sum1
    end
    #Put x,y polynomials together
    for i = 1:4
        for j = 1:4
            #how can I put all these variables together symbolically?
            sumFunc = sumFunc + v11[i,j]*xPoly[i]*yPoly[j]
        end
    end
    return sumFunc
end

function plotty(sumFunc,x,y,v1,ele)
    xx = makeGrid(x,ele)
    yy = makeGrid(y,ele)
    v11 = makeGrid(v1,ele)
    homes = zeros[4,4]
    for i= 1:4
        for j= 1:4
            #how can I plug values into the symbolic Lagrange polynomial?
            homes[i,j] = sumFunc(xx[i,j],yy[i,j])
        end
    end
    plot(xx,yy,homes)
end

#To plot each element of v1 individually
for i = 1:1024
    #is there some way to concatinate graphs so that the functions can be overlayed on each other?
    sumFunc = v1_fit(x,y,v1,i)
    plotty(sumFunc,x,y,v1,i)
end

#To plot each element of v2 individually
for i = 1:1024
    sumFunc = v1_fit(x,y,v2,i)
    plotty(sumFunc,x,y,v2,i)
end


