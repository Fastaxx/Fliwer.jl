using Fliwer
using IterativeSolvers
using WriteVTK

### 2D Test Case : Diphasic Unsteady Diffusion Equation with a Disk
# Define the mesh
nx, ny = 320, 320
lx, ly = 8., 8.
x0, y0 = 0., 0.
domain = ((x0, lx), (y0, ly))
mesh = CartesianMesh((nx, ny), (lx, ly), (x0, y0))

# Define the body
radius, center = ly/4, (lx/2, ly/2)
circle = Body((x,y,_=0)->sqrt((x-center[1])^2 + (y-center[2])^2) - radius, (x,y,_)->(x,y), domain, false)
circle_c = Body((x,y,_=0)->-(sqrt((x-center[1])^2 + (y-center[2])^2) - radius), (x,y,_)->(x,y), domain, false)

# Identify cells
identify!(mesh, circle)

# Define the capacity
capacity = Capacity(circle, mesh)
capacity_c = Capacity(circle_c, mesh)

# Define the operators
operator = DiffusionOps(capacity.A, capacity.B, capacity.V, capacity.W, (nx+1, ny+1))
operator_c = DiffusionOps(capacity_c.A, capacity_c.B, capacity_c.V, capacity_c.W, (nx+1, ny+1))

# Define the boundary conditions
bc = Dirichlet(0.0)
bc_b = BorderConditions(Dict{Symbol, AbstractBoundary}(:left => bc, :right => bc, :top => bc, :bottom => bc))

ic = InterfaceConditions(ScalarJump(1.0, 2.0, 0.0), FluxJump(1.0, 1.0, 0.0))

# Define the source term
f1 = (x,y,z,t)->0.0
f2 = (x,y,z,t)->0.0

# Define the phases
Fluide_1 = Phase(capacity, operator, f1, 1.0)
Fluide_2 = Phase(capacity_c, operator_c, f2, 1.0)

# Initial condition
u0ₒ1 = ones((nx+1)*(ny+1))
u0ᵧ1 = ones((nx+1)*(ny+1))
u0ₒ2 = zeros((nx+1)*(ny+1))
u0ᵧ2 = ones((nx+1)*(ny+1))
u0 = vcat(u0ₒ1, u0ᵧ1, u0ₒ2, u0ᵧ2)

# Define the solver
Δt = 0.01
Tend = 1.0
solver = DiffusionUnsteadyDiph(Fluide_1, Fluide_2, bc_b, ic, Δt, Tend, u0, "BE")

# Solve the problem
Fliwer.solve_DiffusionUnsteadyDiph!(solver, Fluide_1, Fluide_2, u0, Δt, Tend, bc_b, ic, "CN"; method=Base.:\)
#Fliwer.solve_DiffusionUnsteadyDiph!(solver, Fluide_1, Fluide_2, u0, Δt, Tend, bc_b, ic; method=IterativeSolvers.gmres, maxiter=10000, verbose=false)

# Write the solution to a VTK file
#write_vtk("solution", mesh, solver)

# Plot the solution
plot_solution(solver, mesh, circle, capacity, state_i=101)

readline()
# Plot the Profile
#plot_profile(solver, mesh; x=lx/2.01)

# Animation
animate_solution(solver, mesh, circle)

# Analytical solution
using QuadGK
using SpecialFunctions

Dg,Dl = 1.0, 1.0
R0 = radius
cg0, cl0 = 1.0, 0.0
He = 0.5

D = sqrt(Dg/Dl)

function Phi(u)
    term1 = Dg*sqrt(Dl)*besselj1(u*R0)*bessely0(D*u*R0)
    term2 = He*Dl*sqrt(Dg)*besselj0(u*R0)*bessely1(D*u*R0)
    return term1 - term2
end

function Psi(u)
    term1 = Dg*sqrt(Dl)*besselj1(u*R0)*besselj0(D*u*R0)
    term2 = He*Dl*sqrt(Dg)*besselj0(u*R0)*besselj1(D*u*R0)
    return term1 - term2
end

function cg_integrand(u, r, t)
    Φu = Phi(u)
    Ψu = Psi(u)
    denom = u^2*(Φu^2 + Ψu^2)
    num   = exp(-Dg*u^2*t)*besselj0(u*r)*besselj1(u*R0)
    return iszero(denom) ? 0.0 : num/denom
end

function cl_integrand(u, r, t)
    Φu = Phi(u)
    Ψu = Psi(u)
    denom = u*(Φu^2 + Ψu^2)
    term1 = besselj0(D*u*r)*Φu
    term2 = bessely0(D*u*r)*Ψu
    num   = exp(-Dg*u^2*t)*besselj1(u*R0)*(term1 - term2)
    return iszero(denom) ? 0.0 : num/denom
end

function compute_cg(r_values, t_values)
    prefactor = (4*cg0*Dg*Dl*Dl*He)/(π^2*R0)
    cg_results = Matrix{Float64}(undef, length(t_values), length(r_values))
    for (i, t) in pairs(t_values)
        Umax = 5.0/sqrt(Dg*t)
        for (j, r) in pairs(r_values)
            val, _ = quadgk(u->cg_integrand(u, r, t), 0, Umax; atol=1e-6, rtol=1e-6)
            cg_results[i, j] = prefactor*val
        end
    end
    return cg_results
end

function compute_cl(r_values, t_values)
    prefactor = (2*cg0*Dg*sqrt(Dl)*He)/π
    cl_results = Matrix{Float64}(undef, length(t_values), length(r_values))
    for (i, t) in pairs(t_values)
        Umax = 5.0/sqrt(Dg*t)
        for (j, r) in pairs(r_values)
            val, _ = quadgk(u->cl_integrand(u, r, t), 0, Umax; atol=1e-6, rtol=1e-6)
            cl_results[i, j] = prefactor*val
        end
    end
    return cl_results
end

r_values_inside = range(1e-6, stop=R0, length=100)
r_values_outside = range(R0, stop=4*R0, length=100)
r_values = range(1e-6, stop=4*R0, length=nx+1)
t_values = [1.0]

cg_vals = compute_cg(collect(r_values), t_values)
cl_vals = compute_cl(collect(r_values), t_values)

# Plot the analytical solution and the numerical solution (profile)
x=2.01
y = range(mesh.x0[2], stop=mesh.x0[2]+mesh.h[2][1]*length(mesh.h[2]), length=length(mesh.h[2])+1)
x_idx = round(Int, (x - mesh.x0[1]) / mesh.h[1][1])
state = solver.states[end]
u1ₒ = reshape(state[1:length(state) ÷ 4], (length(mesh.centers[1])+1, length(mesh.centers[2])+1))[x_idx, 1:div(length(mesh.centers[2])+1, 2)]
u2ₒ = reshape(state[2*length(state) ÷ 4 + 1:3*length(state) ÷ 4], (length(mesh.centers[1])+1, length(mesh.centers[2])+1))[x_idx, 1:div(length(mesh.centers[2])+1, 2)]

using CairoMakie

fig = Figure()
ax = Axis(fig[1, 1], xlabel="r", ylabel="c", title="Analytical solution")
lines!(ax, r_values, cg_vals[1, :], color=:blue, linewidth=2, label="Analytical solution - Phase 1")
scatter!(ax, r_values, u1ₒ, color=:red, label="Numerical solution - Phase 1")
lines!(ax, r_values, cl_vals[1, :], color=:red, linewidth=2, label="Analytical solution - Phase 2")
scatter!(ax, r_values,u2ₒ, color=:blue, label="Numerical solution - Phase 2")
axislegend(ax)
display(fig)
