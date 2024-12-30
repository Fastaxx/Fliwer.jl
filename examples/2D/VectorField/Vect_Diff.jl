using Fliwer
using IterativeSolvers

### 2D Test Case : Monophasic Unsteady Diffusion of a Vector Field inside a Disk
# Define the mesh
nx, ny = 20, 40
lx, ly = 1., 1.
x0, y0 = 0., 0.
domain = ((x0, lx), (y0, ly))
mesh = CartesianMesh((nx, ny), (lx, ly), (x0, y0))

# two new staggered grids are defined, each corresponding to a velocity component, one displaced half a cell width and the other one displaced half a cell height.
mesh_u = CartesianMesh((nx-1, ny), (lx - lx/nx, ly), (x0 + lx/(2*nx), y0))
mesh_v = CartesianMesh((nx, ny-1), (lx, ly - ly/ny), (x0, y0 + ly/(2*ny)))

@show size(mesh_u.centers[1]), size(mesh_u.centers[2])
@show size(mesh_v.centers[1]), size(mesh_v.centers[2])
# Define the body
radius, center = ly/4, (lx/2, ly/2) #.+ (0.01, 0.01)
circle = Body((x,y,_=0)->(sqrt((x-center[1])^2 + (y-center[2])^2) - radius), (x,y,_)->(x,y), domain, false)

# Identify cells
identify!(mesh, circle)
identify!(mesh_u, circle)
identify!(mesh_v, circle)

# Define the capacity
capacity = Capacity(circle, mesh)
capacity_u = Capacity(circle, mesh_u)
capacity_v = Capacity(circle, mesh_v)

# Define the operators
operator_u = DiffusionOps(capacity_u.A, capacity_u.B, capacity_u.V, capacity_u.W, (nx, ny+1))
operator_v = DiffusionOps(capacity_v.A, capacity_v.B, capacity_v.V, capacity_v.W, (nx+1, ny))

# Define the boundary conditions for each velocity component
# For u-component
bc_u = BorderConditions(Dict{Symbol, AbstractBoundary}(:left => Dirichlet(0.0), :right => Dirichlet(0.0), :top => Dirichlet(0.0), :bottom => Dirichlet(0.0)))
ic_u = Robin(1.0, 1.0, 0.0)
# For v-component
bc_v = BorderConditions(Dict{Symbol, AbstractBoundary}(:left => Dirichlet(0.0), :right => Dirichlet(0.0), :top => Dirichlet(0.0), :bottom => Dirichlet(0.0)))
ic_v = Dirichlet(0.0)

# Define the source term
fu = (x,y,z,t)-> 1.0
fv = (x,y,z,t)-> 1.0

# Define the phase
Fluide = VectorPhase((capacity_u, capacity_v), (operator_u, operator_v), (fu,fv), 1.0)

# Initial condition
uxₒ0 = zeros(nx*(ny+1))
uyₒ0 = zeros((nx+1)*ny)

uxᵧ0 = zeros(nx*(ny+1))
uyᵧ0 = zeros((nx+1)*ny)

u0 = vcat(uxₒ0, uxᵧ0, uyₒ0, uyᵧ0)

u0x = (uxₒ0, uxᵧ0)
u0y = (uyₒ0, uyᵧ0)

# Define the solver
Δt = 0.01
Tend = 1.0
solver = DiffusionVecUnsteadyMono(Fluide, (bc_u, bc_v), (ic_u, ic_v), Δt, Tend, u0x, u0y)

# Solve the problem
solve_DiffusionVecUnsteadyMono!(solver, Fluide, u0x, u0y, Δt, Tend, (bc_u, bc_v), (ic_u, ic_v), method=Base.:\)

write_vtk("DiffusionVecUnsteadyMono", mesh, solver)

# Plot the solution
plot_solution_vector(solver, mesh, circle, capacity_u)