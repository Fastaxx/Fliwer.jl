using Fliwer
using IterativeSolvers

### 3D Test Case : Monophasic Unsteady Diffusion of a Vector Field inside a Sphere
# Define the mesh
nx, ny, nz = 10, 11, 12
lx, ly, lz = 1., 1., 1.
x0, y0, z0 = 0., 0., 0.
domain = ((x0, lx), (y0, ly), (z0, lz))
mesh = CartesianMesh((nx, ny, nz), (lx, ly, lz), (x0, y0, z0))

# three new staggered grids are defined, each corresponding to a velocity component, one displaced half a cell width, the other one displaced half a cell height and the last one displaced half a cell depth.
mesh_u = CartesianMesh((nx-1, ny, nz), (lx - lx/nx, ly, lz), (x0 + lx/(2*nx), y0, z0))
mesh_v = CartesianMesh((nx, ny-1, nz), (lx, ly - ly/ny, lz), (x0, y0 + ly/(2*ny), z0))
mesh_w = CartesianMesh((nx, ny, nz-1), (lx, ly, lz - lz/nz), (x0, y0, z0 + lz/(2*nz)))

# Define the body
radius, center = ly/4, (lx/2, ly/2, lz/2) #.+ (0.01, 0.01)
sphere = Body((x,y,z,_=0)->(sqrt((x-center[1])^2 + (y-center[2])^2 + (z-center[3])^2) - radius), (x,y,z,_)->(x,y,z), domain, false)

# Identify cells
identify!(mesh, sphere)
identify!(mesh_u, sphere)
identify!(mesh_v, sphere)
identify!(mesh_w, sphere)

# Define the capacity
capacity = Capacity(sphere, mesh)
capacity_u = Capacity(sphere, mesh_u)
capacity_v = Capacity(sphere, mesh_v)
capacity_w = Capacity(sphere, mesh_w)

# Define the operators
operator_u = DiffusionOps(capacity_u.A, capacity_u.B, capacity_u.V, capacity_u.W, (nx, ny+1, nz+1))
operator_v = DiffusionOps(capacity_v.A, capacity_v.B, capacity_v.V, capacity_v.W, (nx+1, ny, nz+1))
operator_w = DiffusionOps(capacity_w.A, capacity_w.B, capacity_w.V, capacity_w.W, (nx+1, ny+1, nz))

# Define the boundary conditions for each velocity component
# For u-component
bc_u = BorderConditions(Dict{Symbol, AbstractBoundary}(:left => Dirichlet(0.0), :right => Dirichlet(0.0), :top => Dirichlet(0.0), :bottom => Dirichlet(0.0), :front => Dirichlet(0.0), :back => Dirichlet(0.0)))
ic_u = Robin(1.0, 1.0, 0.0)

# For v-component
bc_v = BorderConditions(Dict{Symbol, AbstractBoundary}(:left => Dirichlet(0.0), :right => Dirichlet(0.0), :top => Dirichlet(0.0), :bottom => Dirichlet(0.0), :front => Dirichlet(0.0), :back => Dirichlet(0.0)))
ic_v = Dirichlet(0.0)

# For w-component
bc_w = BorderConditions(Dict{Symbol, AbstractBoundary}(:left => Dirichlet(0.0), :right => Dirichlet(0.0), :top => Dirichlet(0.0), :bottom => Dirichlet(0.0), :front => Dirichlet(0.0), :back => Dirichlet(0.0)))
ic_w = Dirichlet(0.0)

# Define the source term
fu = (x,y,z,t)-> 1.0
fv = (x,y,z,t)-> 1.0
fw = (x,y,z,t)-> 1.0

# Define the phase
Fluide = VectorPhase((capacity_u, capacity_v, capacity_w), (operator_u, operator_v, operator_w), (fu,fv,fw), 1.0)

# Initial condition
uxₒ0 = zeros(nx*(ny+1)*(nz+1))
uyₒ0 = zeros((nx+1)*ny*(nz+1))
uzₒ0 = zeros((nx+1)*(ny+1)*nz)

uxᵧ0 = zeros(nx*(ny+1)*(nz+1))
uyᵧ0 = zeros((nx+1)*ny*(nz+1))
uzᵧ0 = zeros((nx+1)*(ny+1)*nz)

u0 = vcat(uxₒ0, uxᵧ0, uyₒ0, uyᵧ0, uzₒ0, uzᵧ0)

u0x = (uxₒ0, uxᵧ0)
u0y = (uyₒ0, uyᵧ0)
u0z = (uzₒ0, uzᵧ0)

# Define the solver
Δt = 0.01
Tend = 1.0
solver = DiffusionVecUnsteadyMono(Fluide, (bc_u, bc_v, bc_w), (ic_u, ic_v, ic_w), Δt, Tend, u0x, u0y, u0z)

# Solve the problem
solve_DiffusionVecUnsteadyMono!(solver, Fluide, u0x, u0y, u0z, Δt, Tend, (bc_u, bc_v, bc_w), (ic_u, ic_v, ic_w), method=IterativeSolvers.gmres)
