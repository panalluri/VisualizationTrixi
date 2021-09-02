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

scatter(vec(x), vec(y), vec(vort), zcolor=vec(vort), markerstrokewidth=0, ratio = 1, cam=(0,90))
#scatter(vec(x), vec(y), vec(dudy), zcolor=vec(dudy), markerstrokewidth=0, ratio = 1, cam=(0,90))