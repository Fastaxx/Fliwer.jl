using Fliwer
using IterativeSolvers

### 2D Test Case : Monophasic Unsteady Diffusion Equation inside a moving Disk
# Define the mesh
nx, ny = 80, 80
lx, ly = 4., 4.
x0, y0 = 0., 0.
domain = ((x0, lx), (y0, ly))
mesh = CartesianMesh((nx, ny), (lx, ly), (x0, y0))

# Define the time mesh
Δt = 0.01
Tend = 0.1
nt = Int(Tend/Δt)
t = [i*Δt for i in 0:nt]

# Define the body
radius, center = ly/4, (lx/2, ly/2) .+ (0.01, 0.01)
mapping = (x, y, t) -> (x + 0.1 * t, y)
initial_body = Body((x, y,_=0) -> sqrt((x - center[1])^2 + (y - center[2])^2) - radius, (x,y,_)->(x,y), domain, false)
body = Body((x, y, t) -> sqrt((x - mapping(x, y, t)[1])^2 + (y - mapping(x, y, t)[2])^2) - radius, (x,y,_)->(x,y), domain, false)

# Define the space-time mesh
spaceTimeMesh = CartesianSpaceTimeMesh(mesh, t[1:2])

# Identify cells
identify!(mesh, initial_body)
spaceTimeMesh.tag = mesh.tag

# Define the capacity
capacity = Capacity(body, spaceTimeMesh)
capacity_init = Capacity(initial_body, mesh)

# Define the operators
operator = SpaceTimeOps(capacity.A, capacity.B, capacity.V, capacity.W, (nx+1, ny+1, 2))
operator_init = DiffusionOps(capacity_init.A, capacity_init.B, capacity_init.V, capacity_init.W, (nx+1, ny+1))

# Define the boundary conditions
bc = Dirichlet(0.0)
bc1 = Dirichlet(1.0)

bc_b = BorderConditions(Dict{Symbol, AbstractBoundary}(:left => bc, :right => bc, :top => bc, :bottom => bc))

# Define the source term
f = (x,y,z,t)-> 1.0 #sin(x)*cos(10*y)

Fluide = Phase(capacity, operator, f, 1.0)

# Initial condition
u0ₒ = ones((nx+1)*(ny+1))
u0ᵧ = ones((nx+1)*(ny+1))
u0 = vcat(u0ₒ, u0ᵧ)

# Define the solver
solver = MovingDiffusionUnsteadyMono(Fluide, bc_b, bc1, Δt, Tend, u0, "CN")

# Solve the problem
solve_MovingDiffusionUnsteadyMono!(solver, Fluide, u0, Δt, Tend, nt, bc_b, bc1, body, mesh, t, "CN"; method=Base.:\)

# Plot
using CairoMakie

xₒ = solver.x[1:end÷2]
xᵧ = solver.x[end÷2+1:end]

xₒ = reshape(xₒ, (nx+1, ny+1))
xᵧ = reshape(xᵧ, (nx+1, ny+1))

fig = Figure()
ax = Axis(fig[1, 1])
heatmap!(ax, xₒ, colormap=:viridis)
display(fig)

fig = Figure()
ax = Axis(fig[1, 1])
heatmap!(ax, xᵧ, colormap=:viridis)
display(fig)




# Plot the solution
plot_solution(solver, mesh, body, capacity; state_i=1)

# Animation
#animate_solution(solver, mesh, body)
