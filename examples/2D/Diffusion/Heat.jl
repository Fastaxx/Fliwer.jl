using Fliwer
using IterativeSolvers

### 2D Test Case : Monophasic Unsteady Diffusion Equation inside a Disk
# Define the mesh
nx, ny = 60, 60
lx, ly = 4., 4.
x0, y0 = 0., 0.
domain = ((x0, lx), (y0, ly))
mesh = CartesianMesh((nx, ny), (lx, ly), (x0, y0))

# Define the body
radius, center = ly/4, (lx/2, ly/2) .+ (0.01, 0.01)
circle = Body((x,y,_=0)->(sqrt((x-center[1])^2 + (y-center[2])^2) - radius), (x,y,_)->(x,y), domain, false)

# Identify cells
identify!(mesh, circle)

# Define the capacity
capacity = Capacity(circle, mesh)

# Define the operators
operator = DiffusionOps(capacity.A, capacity.B, capacity.V, capacity.W, (nx+1, ny+1))

cell_types = capacity.cell_types

# Define the boundary conditions 
bc = Dirichlet(1.0)
bc0 = Dirichlet(0.0)

bc_b = BorderConditions(Dict{Symbol, AbstractBoundary}(:left => bc0, :right => bc0, :top => bc0, :bottom => bc0))

# Define the source term
f = (x,y,z,t)->0.0

# Define the phase
Fluide = Phase(capacity, operator, f, 1.0)

# Initial condition
u0ₒ = zeros((nx+1)*(ny+1))
u0ᵧ = zeros((nx+1)*(ny+1))
u0 = vcat(u0ₒ, u0ᵧ)

# Define the solver
Δt = 0.01
Tend = 1.0
solver = DiffusionUnsteadyMono(Fluide, bc_b, bc, Δt, Tend, u0, "BE") # Start by a backward Euler scheme to prevent oscillation due to CN scheme

# Solve the problem
Fliwer.solve_DiffusionUnsteadyMono!(solver, Fluide, u0, Δt, Tend, bc_b, bc, "CN"; method=Base.:\)
#Fliwer.solve_DiffusionUnsteadyMono!(solver, Fluide, u0, Δt, Tend, bc_b, bc; method=IterativeSolvers.bicgstabl, reltol=1e-40, verbose=false)

# Write the solution to a VTK file
#write_vtk("heat", mesh, solver)

# Plot the solution
plot_solution(solver, mesh, circle, capacity; state_i=10)

# Animation
animate_solution(solver, mesh, circle)
