using Fliwer
using IterativeSolvers

### 1D Test Case : Diphasic Unsteady Diffusion Equation 
# Define the mesh
nx = 80
lx = 8.0
x0 = 0.0
domain=((x0,lx),)
mesh = CartesianMesh((nx,), (lx,), (x0,))

# Define the body
xint = 4.0 + 0.1
body = Body((x, _=0) -> (x - xint),(x,_=0)->(x),domain,false)
body_c = Body((x, _=0) -> -(x - xint),(x,_=0)->(x),domain,false)

# Identify cells
identify!(mesh, body)

# Define the capacity
capacity = Capacity(body, mesh)
capacity_c = Capacity(body_c, mesh)

# Define the operators
operator = DiffusionOps(capacity.A, capacity.B, capacity.V, capacity.W, (nx+1,))
operator_c = DiffusionOps(capacity_c.A, capacity_c.B, capacity_c.V, capacity_c.W, (nx+1,))

# Define the boundary conditions
bc1 = Dirichlet(0.0)
bc0 = Dirichlet(1.0)
bc_b = BorderConditions(Dict{Symbol, AbstractBoundary}(:top => bc0, :bottom => bc1))

ic = InterfaceConditions(ScalarJump(1.0, 0.75, 0.0), FluxJump(1.0, 4.0, 0.0))

# Define the source term
f1 = (x,y,z,t)->0.0
f2 = (x,y,z,t)->0.0

# Define the phases
Fluide_1 = Phase(capacity, operator, f1, 1.0)
Fluide_2 = Phase(capacity_c, operator_c, f2, 1.0)

# Initial condition
u0ₒ1 = zeros(nx+1)
u0ᵧ1 = zeros(nx+1)
u0ₒ2 = ones(nx+1)
u0ᵧ2 = ones(nx+1)

u0 = vcat(u0ₒ1, u0ᵧ1, u0ₒ2, u0ᵧ2)

# Define the solver
Δt = 0.01
Tend = 1.0
solver = DiffusionUnsteadyDiph(Fluide_1, Fluide_2, bc_b, ic, Δt, Tend, u0)

# Solve the problem
solve_DiffusionUnsteadyDiph!(solver, Fluide_1, Fluide_2, u0, Δt, Tend, bc_b, ic; method=Base.:\)

# Write the solution to a VTK file
#write_vtk("solution", mesh, solver)

# Plot the solution
plot_solution(solver, mesh, body, capacity)

# Animation
#animate_solution(solver, mesh, body)

# Plot the solution
state_i = 10
u1ₒ = solver.states[state_i][1:nx+1]
u1ᵧ = solver.states[state_i][nx+2:2*(nx+1)]
u2ₒ = solver.states[state_i][2*(nx+1)+1:3*(nx+1)]
u2ᵧ = solver.states[state_i][3*(nx+1)+1:end]

x = range(x0, stop = lx, length = nx+1)
using CairoMakie

fig = Figure()
ax = Axis(fig[1, 1], xlabel="x", ylabel="u", title="Diphasic Unsteady Diffusion Equation")

for (idx, i) in enumerate(1:10:length(solver.states))
    u1ₒ = solver.states[i][1:nx+1]
    u1ᵧ = solver.states[i][nx+2:2*(nx+1)]
    u2ₒ = solver.states[i][2*(nx+1)+1:3*(nx+1)]
    u2ᵧ = solver.states[i][3*(nx+1)+1:end]

    # Display labels only for first iteration
    scatter!(ax, x, u1ₒ, color=:blue,  markersize=3,
        label = (idx == 1 ? "Bulk Field - Phase 1" : nothing))
    scatter!(ax, x, u1ᵧ, color=:red,   markersize=3,
        label = (idx == 1 ? "Interface Field - Phase 1" : nothing))
    scatter!(ax, x, u2ₒ, color=:green, markersize=3,
        label = (idx == 1 ? "Bulk Field - Phase 2" : nothing))
    scatter!(ax, x, u2ᵧ, color=:orange, markersize=3,
        label = (idx == 1 ? "Interface Field - Phase 2" : nothing))
end

axislegend(ax, position=:rb)
display(fig)