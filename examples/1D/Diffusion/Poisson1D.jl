using Fliwer
using IterativeSolvers
using LinearAlgebra

### 1D Test Case : Monophasic Steady Diffusion Equation
# Define the mesh
nx = 160
lx = 1.
x0 = 0.
domain = ((x0, lx),)
mesh = CartesianMesh((nx,), (lx,), (x0,))

# Define the body
pos = 0.5
body = Body((x,_=0)->(x - pos), (x,)->(x,), domain, false)

# Identify cells
identify!(mesh, body)

# Define the capacity
capacity = Capacity(body, mesh)

# Define the operators
operator = DiffusionOps(capacity.A, capacity.B, capacity.V, capacity.W, (nx+1,))

# Define the boundary conditions
bc = Dirichlet(1.0)
bc1 = Dirichlet(0.0)

bc_b = BorderConditions(Dict{Symbol, AbstractBoundary}(:top => bc1, :bottom => bc1))

# Define the source term
f = (x, y, _=0) -> 0.0

Fluide = Phase(capacity, operator, f, 1.0)

# Define the solver
solver = DiffusionSteadyMono(Fluide, bc_b, bc)

# Solve the problem
solve_DiffusionSteadyMono!(solver, Fluide; method=Base.:\)

# Plot the solution
plot_solution(solver, mesh, body, capacity)

# Analytical solution
u_analytical = (x) -> (x < pos) ? 2x : 0.0

err, global_err, full_err, cut_err, empty_err = check_convergence((x) -> (x < pos) ? 2x : 0.0, solver, capacity, 2)

# Plot the error
using CairoMakie

fig = Figure(size = (800, 600))
ax = Axis(fig[1, 1], xlabel = "x", ylabel = "u(x) - u_num(x)")
scatter!(ax, err, color = :blue)
fig