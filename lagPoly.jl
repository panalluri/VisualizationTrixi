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

xVec = x[1,:]
append!(x[1,:],[1])
yVec=xVec

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
end
