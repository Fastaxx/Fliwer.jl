using Fliwer
using IterativeSolvers

### 2D Test Case : Diphasic Steady Advection-Diffusion Equation inside a Disk
# Define the mesh
nx, ny = 40, 40
lx, ly = 4., 4.
x0, y0 = 0., 0.
domain = ((x0, lx), (y0, ly))
mesh = CartesianMesh((nx, ny), (lx, ly), (x0, y0))

# Define the body
radius, center = ly/4, (lx/2, ly/2) .+ (0.1, 0.1)
circle = Body((x,y,_=0)->sqrt((x-center[1])^2 + (y-center[2])^2) - radius, (x,y,_)->(x,y), domain, false)
circle_c = Body((x,y,_=0)->-(sqrt((x-center[1])^2 + (y-center[2])^2) - radius), (x,y,_)->(x,y), domain, false)

# Identify cells
identify!(mesh, circle)

# Define the capacity
capacity = Capacity(circle, mesh)
capacity_c = Capacity(circle_c, mesh)

# Initialize the velocity field by setting the velocity to zero
uₒx, uₒy = zeros((nx+1)*(ny+1)), zeros((nx+1)*(ny+1))
uγx, uγy = zeros((nx+1)*(ny+1)), zeros((nx+1)*(ny+1))

uₒ, uᵧ = vcat(uₒx, uₒy), vcat(uγx, uγy)

# Define the operators
operator = ConvectionOps(capacity.A, capacity.B, capacity.V, capacity.W, (nx+1, ny+1), uₒ, uᵧ)
operator_c = ConvectionOps(capacity_c.A, capacity_c.B, capacity_c.V, capacity_c.W, (nx+1, ny+1), uₒ, uᵧ)

# Define the boundary conditions
bc = Dirichlet(0.0)
bc_b = BorderConditions(Dict{Symbol, AbstractBoundary}(:left => bc, :right => bc, :top => bc, :bottom => bc))

ic = InterfaceConditions(ScalarJump(1.0, 1.0, 0.0), FluxJump(1.0, 1.0, 0.0))

# Define the source term
f1 = (x,y,_)->1.0 #cos(x)*sin(10*y)
f2 = (x,y,_)->0.0 #cos(x)*sin(10*y)

# Define the phases
Fluide_1 = Phase(capacity, operator, f1, 1.0)
Fluide_2 = Phase(capacity_c, operator_c, f2, 1.0)

# Define the solver
solver = AdvectionDiffusionSteadyDiph(Fluide_1, Fluide_2, bc_b, ic)

# Solve the problem
solve_AdvectionDiffusionSteadyDiph!(solver, Fluide_1, Fluide_2; method=IterativeSolvers.gmres, verbose=false, reltol=1e-40)

# Plot the solution usign Makie
plot_solution(solver, mesh, circle, capacity)