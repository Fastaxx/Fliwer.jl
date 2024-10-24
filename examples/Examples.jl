using Fliwer

### 2D Test Case : Monophasic Steady Diffusion Equation inside a Disk
# Define the mesh
nx, ny = 40, 40
lx, ly = 4., 4.
x0, y0 = 0., 0.
domain = ((x0, lx), (y0, ly))
mesh = CartesianMesh((nx, ny), (lx, ly), (x0, y0))

# Define the body
radius, center = ly/4, (lx/2, ly/2)
circle = Body((x,y,_=0)->sqrt((x-center[1])^2 + (y-center[2])^2) - radius, (x,y,_)->(x,y), domain, false)

# Define the capacity
capacity = Capacity(circle, mesh)

# Define the operators
operator = DiffusionOps(capacity.A, capacity.B, capacity.V, capacity.W, (nx, ny))

# Define the boundary conditions
bc = Dirichlet(0.0)

# Define the source term
f = (x,y,_)->1.0

# Define the solver
solver = DiffusionSolver(true, true, 1.0)

# Initialize the solver
initialize!(solver,operator,bc)

# Solve the problem
u = solve!(solver, operator, bc, f)

# Plot the solution

using Plots
Plots.default(show=true)
heatmap(reshape(u[1:length(u)÷2],(nx,ny))',colormap=:viridis)
readline()

println(maximum(u))