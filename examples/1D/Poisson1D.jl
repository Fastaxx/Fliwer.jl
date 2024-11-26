using Fliwer
using IterativeSolvers

### 1D Test Case : Monophasic Steady Diffusion Equation
# Define the mesh
nx = 20
lx = 4.
x0 = 0.
domain = ((x0, lx),)
mesh = CartesianMesh((nx,), (lx,), (x0,))

# Define the body
pos = 2. + 0.1
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
solve!(solver, Fluide; method=IterativeSolvers.bicgstabl, verbose=false)

# Plot the solution
plot_solution(solver, mesh, body, capacity)