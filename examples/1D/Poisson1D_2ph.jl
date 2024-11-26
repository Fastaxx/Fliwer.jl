using Fliwer
using IterativeSolvers

### 1D Test Case : Diphasic Steady Diffusion Equation
# Define the mesh
nx = 20
lx = 4.0
x0 = 0.
domain = ((x0, lx),)
mesh = CartesianMesh((nx,), (lx,), (x0,))

# Define the body
pos = 2.0+0.1
body = Body((x,_=0)->(x - pos), (x,)->(x,), domain, false)
body_c = Body((x,_=0)->-(x - pos), (x,)->(x,), domain, false)

# Identify cells
identify!(mesh, body)

# Define the capacity
capacity = Capacity(body, mesh)
capacity_c = Capacity(body_c, mesh)

# Define the operators
operator = DiffusionOps(capacity.A, capacity.B, capacity.V, capacity.W, (nx+1,))
operator_c = DiffusionOps(capacity_c.A, capacity_c.B, capacity_c.V, capacity_c.W, (nx+1,))

# Define the boundary conditions
bc = Dirichlet(1.0)
bc1 = Dirichlet(0.0)

bc_b = BorderConditions(Dict{Symbol, AbstractBoundary}(:top => bc, :bottom => bc1))

ic = InterfaceConditions(ScalarJump(1.0, 1.0, 0.0), FluxJump(1.0, 1.0, 0.0))

# Define the source term
f1 = (x, y, _=0) -> 0.0
f2 = (x, y, _=0) -> 0.0

# Define the phases
Fluide_1 = Phase(capacity, operator, f1, 1.0)
Fluide_2 = Phase(capacity_c, operator_c, f2, 1.0)

# Define the solver 
solver = DiffusionSteadyDiph(Fluide_1, Fluide_2, bc_b, ic)

# Solve the problem
solve!(solver, Fluide_1, Fluide_2; method=IterativeSolvers.gmres, verbose=false)

# Plot the solution
plot_solution(solver, mesh, body, capacity)
