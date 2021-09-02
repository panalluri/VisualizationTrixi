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

function Li_xBasis(x,xGrid,xi)
    xPolBase = 0
    for j = 1:16
        #how do I use symbolic x?
        if xi != xGrid[j]
            xPolBase = xPolBase + (x-xGrid[j])/(xi-xGrid[j])
        end
    end
    return xPolBase
end

function Lj_yBasis(y,yGrid,yi)
    yPolBase = 0
    for j = 1:16
        #how do I use symbolic y?
        if yi != yGrid[j]
            yPolBase = yPolBase + (y-yGrid[j])/(yi-yGrid[j])
        end
    end
    return yPolBase
end

function Pxy(xGrid,yGrid,v1Grid)
    poly = zeros(16)
    for i = 1:16
        poly1=0
        for j = 1:16
            xPolBase = Li_xBasis(x[i],xGrid,xGrid[j])
            yPolBase = Lj_yBasis(y[i],yGrid,yGrid[j])
            poly1 = poly1 + v1Grid[j]*xPolBase*yPolBase
        end
        poly[i] = poly1
    end
    return poly
end

function plotty(x,y,v1)
    #for i = 1:1
    i=1
        xGrid = x[:,i]
        yGrid = y[:,i]
        v1Grid = v1[:,i]
        poly = Pxy(xGrid,yGrid,v1Grid)
        xx = makeGrid(x,i)
        yy = makeGrid(y,i)
        poly = makeGrid(poly,1)
        plot(xx,yy,poly)
        # plot(xGrid,yGrid,poly)
    #end
end

plotty(x,y,v1)
