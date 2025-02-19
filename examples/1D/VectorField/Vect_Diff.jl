using Fliwer
using IterativeSolvers

### 1D Test Case : Monophasic Unsteady Diffusion of a Vector Field inside a Disk
# Define the mesh
nx = 20
lx = 1.
x0 = 0.
domain = ((x0, lx),)
mesh = CartesianMesh((nx,), (lx,), (x0,))

mesh_u = CartesianMesh((nx-1,), (lx - lx/nx,), (x0 + lx/(2*nx),))

@show size(mesh_u.centers[1])

# Define the body
radius, center = 0.25, 0.5
circle = Body((x,_=0)->(sqrt((x-center)^2) - radius), (x,_)->(x,), domain, false)

# Identify cells
identify!(mesh, circle)
identify!(mesh_u, circle)

# Define the capacity
capacity = Capacity(circle, mesh)
capacity_u = Capacity(circle, mesh_u)

# Define the operators
operator_u = DiffusionOps(capacity_u.A, capacity_u.B, capacity_u.V, capacity_u.W, (nx,))

# Define the boundary conditions for each velocity component
# For u-component
bc_u = BorderConditions(Dict{Symbol, AbstractBoundary}(:left => Dirichlet(0.0), :right => Dirichlet(0.0)))
ic_u = Robin(1.0, 1.0, 0.0)

# Define the source term
fu = (x,y,z,t)-> 1.0

# Define the phase
Fluide = VectorPhase((capacity_u,), (operator_u,), (fu,), 1.0)

# Initial condition
uxₒ0 = zeros(nx)
uxᵧ0 = zeros(nx)

u0 = vcat(uxₒ0, uxᵧ0)

u0x = (uxₒ0, uxᵧ0)

# Solve the problem
Δt = 0.01
Tend = 1.0
solver = DiffusionVecUnsteadyMono(Fluide, (bc_u,), (ic_u,), Δt, Tend, u0x)

# Solve the problem
solve_DiffusionVecUnsteadyMono!(solver, Fluide, u0x, Δt, Tend, (bc_u,), (ic_u,), method=Base.:\)

write_vtk("DiffusionVecUnsteadyMono", mesh, solver)

# Plot the solution
plot_solution_vector(solver, mesh, circle, capacity_u)