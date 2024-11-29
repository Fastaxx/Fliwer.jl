using Fliwer
using IterativeSolvers

### 2D Test Case : Diphasic Unsteady Advection-Diffusion Equation inside a Disk
# Define the mesh
nx, ny = 40, 40
lx, ly = 8., 8.
x0, y0 = 0., 0.
domain = ((x0, lx), (y0, ly))
mesh = CartesianMesh((nx, ny), (lx, ly), (x0, y0))

# Define the body
radius, center = ly/2, (lx/2, ly/2) #.+ (0.1, 0.1)
circle = Body((x,y,_=0)->sqrt((x-center[1])^2 + (y-center[2])^2) - radius, (x,y,_)->(x,y), domain, false)
circle_c = Body((x,y,_=0)->-(sqrt((x-center[1])^2 + (y-center[2])^2) - radius), (x,y,_)->(x,y), domain, false)

# Identify cells
identify!(mesh, circle)

# Define the capacity
capacity = Capacity(circle, mesh)
capacity_c = Capacity(circle_c, mesh)

# Initialize the velocity field by setting the velocity to zero
uₒx, uₒy = ones((nx+1)*(ny+1)), ones((nx+1)*(ny+1))
uγx, uγy = ones((nx+1)*(ny+1)), ones((nx+1)*(ny+1))

uₒ, uᵧ = vcat(uₒx, uₒy), vcat(uγx, uγy)

# Define the operators
operator = ConvectionOps(capacity.A, capacity.B, capacity.V, capacity.W, (nx+1, ny+1), uₒ, uᵧ)
operator_c = ConvectionOps(capacity_c.A, capacity_c.B, capacity_c.V, capacity_c.W, (nx+1, ny+1), uₒ, uᵧ)

# Define the boundary conditions
bc = Dirichlet(0.0)
bc_b = BorderConditions(Dict{Symbol, AbstractBoundary}(:left => bc, :right => bc, :top => bc, :bottom => bc))

ic = InterfaceConditions(ScalarJump(1.0, 1.0, 0.0), FluxJump(1.0, 1.0, 0.0))

# Define the source term
f1 = (x,y,z,t)->0.0
f2 = (x,y,z,t)->0.0

# Define the phases
Fluide_1 = Phase(capacity, operator, f1, 1.0)
Fluide_2 = Phase(capacity_c, operator_c, f2, 1.0)

# Initial condition
u0ₒ1 = ones((nx+1)*(ny+1))
u0ᵧ1 = ones((nx+1)*(ny+1))
u0ₒ2 = zeros((nx+1)*(ny+1))
u0ᵧ2 = zeros((nx+1)*(ny+1))
u0 = vcat(u0ₒ1, u0ᵧ1, u0ₒ2, u0ᵧ2)

# Define the solver
Δt = 0.01
Tend = 1.0
solver = AdvectionDiffusionUnsteadyDiph(Fluide_1, Fluide_2, bc_b, ic, Δt, Tend, u0)

# Solve the problem
solve_AdvectionDiffusionUnsteadyDiph!(solver, Fluide_1, Fluide_2, u0, Δt, Tend, bc_b, ic; method=IterativeSolvers.gmres, maxiter=10000, verbose=false)

# Write the solution to a VTK file
write_vtk("solution", mesh, solver)

# Plot the solution
plot_solution(solver, mesh, circle, capacity)