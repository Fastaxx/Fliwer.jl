using Fliwer
using IterativeSolvers

### 2D Test Case : Diphasic Steady Diffusion Equation inside a Disk
# Define the mesh
nx, ny = 80, 80
lx, ly = 4., 4.
x0, y0 = 0., 0.
domain = ((x0, lx), (y0, ly))
mesh = CartesianMesh((nx, ny), (lx, ly), (x0, y0))

# Define the body
radius, center = ly/4, (lx/2, ly/2) .+ (0.01, 0.01)
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
bc = Dirichlet(1.0)
bc1 = Dirichlet(0.0)
bc2 = Dirichlet(2.0)
bc_b = BorderConditions(Dict{Symbol, AbstractBoundary}(:left => bc, :right => bc1, :top => bc, :bottom => bc2))

ic = InterfaceConditions(ScalarJump(1.0, 1.0, 0.0), FluxJump(1.0, 1.0, 0.0))

# Define the source term
f1 = (x,y,_)->0.0 #cos(x)*sin(10*y)
f2 = (x,y,_)->0.0 #cos(x)*sin(10*y)

# Define the phases
Fluide_1 = Phase(capacity, operator, f1, 1.0)
Fluide_2 = Phase(capacity_c, operator_c, f2, 1.0)

# Define the solver
solver = DiffusionSteadyDiph(Fluide_1, Fluide_2, bc_b, ic)

# Solve the problem
solve!(solver, Fluide_1, Fluide_2; method=IterativeSolvers.gmres, verbose=false)

# Plot the solution usign Makie
plot_solution(solver, mesh, circle)

# Write the solution to a VTK file
write_vtk("poisson_2d", mesh, solver)