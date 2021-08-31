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
    xPoly = zeros(4)
    yPoly = zeros(4)
    #Find x polynomials
    for i = 1:4
        sum1 = 0
        for j = 1:4
            sum1 = sum1 + (x-xLine[j])/(xLine[i]-xLine[j])
        end
        xPoly[i] = sum1
    end
    #Find y polynomials
    for i = 1:4
        sum1 = 0
        for j = 1:4
            sum1 = sum1 + (y-yLine[j])/(yLine[i]-yLine[j])
        end
        yPoly[i] = sum1
    end
    #Put x,y polynomials together
    for i = 1:4
        for j = 1:4
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
            homes[i,j] = sumFunc(xx[i,j],yy[i,j])
        end
    end
    plot(xx,yy,homes)
end


