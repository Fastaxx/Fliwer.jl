using Fliwer
using IterativeSolvers

### 2D Test Case : Poiseuille Flow
# Define the mesh
nx, ny = 100, 100
lx, ly = 1., 1.
x0, y0 = 0., 0.
domain = ((x0, lx), (y0, ly))
mesh_p = CartesianMesh((nx, ny), (lx, ly), (x0, y0))

# For u-velocity (faces in x-direction)
mesh_u = CartesianMesh((nx+1, ny), (lx, ly), (x0 - lx/(2*nx), y0))
# For v-velocity (faces in y-direction)
mesh_v = CartesianMesh((nx, ny+1), (lx, ly), (x0, y0 - ly/(2*ny)))

# Define the body
body = NoBody2D(domain)

# Identify cells in the meshes
identify!(mesh_p, body)
identify!(mesh_u, body)
identify!(mesh_v, body)

# Define capacities using the new meshes : Aₐ, Bₐ, Vₐ, Wₐ for a ∈ {p, u, v}
capacity_p = Capacity(body, mesh_p)
capacity_u = Capacity(body, mesh_u)
capacity_v = Capacity(body, mesh_v)

# Initialize the velocity field with a rotating field
uₒx, uₒy = zeros((nx+2)*(ny+1)), zeros((nx+1)*(ny+2))
uₒ = (uₒx, uₒy)

# For boundary velocities, if they are zero:
uᵧ = zeros((nx + 2) * (ny + 1) + (nx + 1) * (ny + 2))

# Define the operators
operator = NavierStokesOps(capacity_p, capacity_u, capacity_v, (nx+1, ny+1), uₒ, uᵧ)

# Define the boundary conditions
bc = BorderConditions(Dict{Symbol, AbstractBoundary}(:left => Dirichlet(0.0), :right => Dirichlet(0.0), :top => Dirichlet(0.0), :bottom => Dirichlet(0.0)))

# Define the source term
f = (x, y, z, t) -> 0.0

# Define the phase
ρ = 1.0
Re = 100.0
pressure = Pressure(capacity_p, operator, f, 1.0)
velocity = Velocity(capacity_u, capacity_v, operator, f, ρ, Re)

# Define the solver
Δt = 0.01
Tend = 1.0
solver = NavierStokesUnsteadyMono(pressure, velocity, bc, Δt, Tend, uₒ)

# Solve the problem
solve_NavierStokesUnsteadyMono!(solver, pressure, velocity, uₒ, Δt, Tend, bc; method=IterativeSolvers.bicgstabl, abstol=1e-15, verbose=false)

# Write the solution to a VTK file
