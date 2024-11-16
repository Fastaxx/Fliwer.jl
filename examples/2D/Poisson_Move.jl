using Fliwer
using IterativeSolvers

### 2D Test Case : Monophasic Steady Diffusion Equation inside a moving Disk
# Define the mesh
nx, ny = 80, 80
lx, ly = 4., 4.
x0, y0 = 0., 0.
domain = ((x0, lx), (y0, ly))
mesh = CartesianMesh((nx, ny), (lx, ly), (x0, y0))

# Define the body
radius, center = ly/4, (lx/2, ly/2) .+ (0.01, 0.01)
mapping = (t, x, y) -> (x + 0.1 * t, y)
circle = Body((x,y,t)->-(sqrt((x-center[1])^2 + (y-center[2])^2) - radius), mapping, domain, true)

# Identify cells
identify!(mesh, circle)

# Define the capacity
capacity = Capacity(circle, mesh)

# Define the operators
operator = DiffusionOps(capacity.A, capacity.B, capacity.V, capacity.W, (nx+1, ny+1))

# Define the boundary conditions
bc = Dirichlet(1.0)
bc1 = Dirichlet(0.0)

bc_b = BorderConditions(Dict{Symbol, AbstractBoundary}(:left => bc, :right => bc, :top => bc1, :bottom => bc))

# Define the source term
f = (x,y,t)-> 4.0 #sin(x)*cos(10*y)

Fluide = Phase(capacity, operator, f, 1.0)

# Define the solver
solver = DiffusionSteadyMono(Fluide, bc_b, bc)

# Solve the problem
solve!(solver, Fluide; method=IterativeSolvers.bicgstabl, verbose=false)

# Plot the solution
plot_solution(solver, mesh, circle)

