using Fliwer

### 2D Test Case : Monophasic Steady Diffusion Equation inside a Disk
# Define the mesh
nx, ny = 20, 20
hx, hy = ones(nx), ones(ny)
x0, y0 = 0., 0.
mesh = CartesianMesh((hx,hy), (x0,y0))

# Define the body
domain = ((0., 1.), (0., 1.))
circle = Body((x,y,_=0)->sqrt(x^2 + y^2) - 10, (x,y,_)->(x,y), domain, false)

# Define the capacity
capacity = Capacity(circle, mesh)

# Define the operators
operator = DiffusionOps(capacity.A, capacity.B, capacity.V, capacity.W, (nx+1, ny+1))

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
heatmap(reshape(u[1:length(u)÷2],(nx+1,ny+1))',colormap=:viridis)
readline()