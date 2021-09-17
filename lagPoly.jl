using Trixi
using StructArrays
using GLMakie

include("C:\\Users\\Prani\\.julia\\dev\\Trixi\\examples\\dgmulti_2d\\elixir_euler_kelvin_helmholtz_instability.jl")
rho, rho_v1, rho_v2, E = StructArrays.components(sol.u[end])
v1 = rho_v1./rho
v2 = rho_v2./rho
@unpack x, y = mesh.md;
x = mesh.md.x
y = mesh.md.y

function compute_derivatives(u, mesh, dg)
    rd = dg.basis # a "RefElemData" rd object
    @unpack rxJ, sxJ, ryJ, syJ, J = mesh.md # unpack out of the "MeshData" md object
    dudx = (rxJ .* (rd.Dr * u) + sxJ .* (rd.Ds * u)) ./ J
    dudy = (ryJ .* (rd.Dr * u) + syJ .* (rd.Ds * u)) ./ J
    return dudx, dudy
end

semi = sol.prob.p
dg = semi.solver
mesh = semi.mesh
dv1dx, dv1dy = compute_derivatives(v1, mesh, dg)
dv2dx, dv2dy = compute_derivatives(v2, mesh, dg)

Vp = dg.basis.Vp
dv2dx_plot = Vp * dv2dx
dv1dy_plot = Vp * dv1dy

vort = dv2dx - dv1dy

# xVec = x[1,:]
# xVec = append!(x[1,:],[1])
# yVec=xVec

#Plot one element

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

function usoln(x,y,xi,yi,xVec,yVec,uij)
    lxi = Li_xBasis(x,xVec,xi)
    lyi = Li_yBasis(y,yVec,yi)
    u_soln = uij*lxi*lyi
    return u_soln
end

function usum(x,y,xVec,yVec,uGrid)
    sum = 0
    for i = 1:length(xVec)
        for j=1:length(yVec)
            lyi = Li_yBasis(y,yVec,y[j])
            lxi = Li_xBasis(x,xVec,x[i])
            sum = sum + uGrid[i,j]*lxi*lyi
        end
    end
    return sum
end

function makeGrid(uVec)
    uGrid = zeros(4,4)
    z = [4 3 2 1]
    for i = 1:4
        uGrid[z[i],:] = transpose(uVec[(i-1)*4+1:i*4])
    end
    return uGrid
end

function elementz(x,y,u,ele)
    xx = x[:,ele]
    yy = sort(y[:,ele])
    xVec = xx[1:4]
    yVec = [yy[1] yy[5] yy[9] yy[13]]
    uVec = u[:,ele]
    return xVec, yVec, uVec
end

function unitPlot(xVec,yVec,uGrid)
    plotty = zeros(length(xVec),length(xVec))
    for i = 1:length(xVec)
        for j = 1:length(yVec)
            plotty[i,j] = usum(x,y,xVec,yVec,uGrid) 
        end
    end
    return plotty
end

ele = 1
xVec, yVec, uVec = elementz(x,y,vort,ele)
uGrid = makeGrid(uVec)
xGrid = makeGrid(x[:,ele])
yGrid = makeGrid(y[:,ele])
plotty = unitPlot(xVec,yVec,uGrid)

plot(xGrid,yGrid,plotty)