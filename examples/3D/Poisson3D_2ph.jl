using Fliwer
using IterativeSolvers

### 3D Test Case : Diphasic Steady Diffusion Equation inside a Sphere
# Define the mesh
nx, ny, nz = 20, 20, 20
lx, ly, lz = 4., 4., 4.
x0, y0, z0 = 0., 0., 0.
domain = ((x0, lx), (y0, ly), (z0, lz))
mesh = CartesianMesh((nx, ny, nz), (lx, ly, lz), (x0, y0, z0))

# Define the body
radius, center = ly/4, (lx/2, ly/2, lz/2) #.+ (0.01, 0.01, 0.01)
sphere = Body((x,y,z)->(sqrt((x-center[1])^2 + (y-center[2])^2 + (z-center[3])^2) - radius), (x,y,z)->(x,y,z), domain, false)
sphere_c = Body((x,y,z)->-(sqrt((x-center[1])^2 + (y-center[2])^2 + (z-center[3])^2) - radius), (x,y,z)->(x,y,z), domain, false)

# Identify cells
identify!(mesh, sphere)

# Define the capacity
capacity = Capacity(sphere, mesh)
capacity_c = Capacity(sphere_c, mesh)

# Define the operators
operator = DiffusionOps(capacity.A, capacity.B, capacity.V, capacity.W, (nx+1, ny+1, nz+1))
operator_c = DiffusionOps(capacity_c.A, capacity_c.B, capacity_c.V, capacity_c.W, (nx+1, ny+1, nz+1))

# Define the boundary conditions
bc = Dirichlet(0.0)
bc_i = Dirichlet(0.0)

bc_b = BorderConditions(Dict{Symbol, AbstractBoundary}(:left => bc, :right => bc, :top => bc, :bottom => bc, :front => bc, :back => bc))

ic = InterfaceConditions(ScalarJump(1.0, 0.5, 0.0), FluxJump(1.0, 1.0, 0.0))

# Define the source term
f1 = (x,y,z)->1.0
f2 = (x,y,z)->1.0

# Define the phases
Fluide_1 = Phase(capacity, operator, f1, 1.0)
Fluide_2 = Phase(capacity_c, operator_c, f2, 1.0)

# Define the solver
solver = DiffusionSteadyDiph(Fluide_1, Fluide_2, bc_b, ic)

# Solve the problem
solve!(solver, Fluide_1, Fluide_2; method=IterativeSolvers.gmres, verbose=false)

# Plot the solution usign Makie
plot_solution(solver, mesh, sphere)