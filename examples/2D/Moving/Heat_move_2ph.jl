using Fliwer
using IterativeSolvers
using LinearAlgebra
using SparseArrays

### 2D Test Case : Diphasic Unsteady Diffusion Equation inside a moving body
# Define the mesh
nx, ny = 20, 20
lx, ly = 8., 8.
x0, y0 = 0., 0.
domain = ((x0, lx), (y0, ly))
mesh = CartesianMesh((nx, ny), (lx, ly), (x0, y0))

# Define the time mesh
Δt = 0.01
Tend = 1.0
nt = Int(round(Tend/Δt))
t = [i*Δt for i in 0:nt]

# Define the body : Translating disk : radius = 0.5, center = (lx/2 + 0.1*t, ly/2)
radius, center = ly/4, (lx/2, ly/2)
c = 0.0

initial_body = Body((x, y,_=0) -> sqrt((x - center[1])^2 + (y - center[2])^2) - radius, (x,y,_)->(x,y), domain, false)
initial_body_c = Body((x, y,_=0) -> -(sqrt((x - center[1])^2 + (y - center[2])^2) - radius), (x,y,_)->(x,y), domain, false)

body = Body((x, y, t) -> sqrt((x - center[1] - c*t)^2 + (y - center[2])^2) - radius, (x,y,_)->(x,y), domain, false)
body_c = Body((x, y, t) -> -(sqrt((x - center[1] - c*t)^2 + (y - center[2])^2) - radius), (x,y,_)->(x,y), domain, false)

# Define the space-time mesh
spaceTimeMesh = CartesianSpaceTimeMesh(mesh, t[1:2])

# Identify cells
identify!(mesh, initial_body)
spaceTimeMesh.tag = mesh.tag

# Define the capacity
capacity = Capacity(body, spaceTimeMesh)
capacity_c = Capacity(body_c, spaceTimeMesh)

capacity_init = Capacity(initial_body, mesh)
capacity_init_c = Capacity(initial_body_c, mesh)

# Define the operators
operator = SpaceTimeOps(capacity.A, capacity.B, capacity.V, capacity.W, (nx+1, ny+1, 2))
operator_c = SpaceTimeOps(capacity_c.A, capacity_c.B, capacity_c.V, capacity_c.W, (nx+1, ny+1, 2))

operator_init = DiffusionOps(capacity_init.A, capacity_init.B, capacity_init.V, capacity_init.W, (nx+1, ny+1))
operator_init_c = DiffusionOps(capacity_init_c.A, capacity_init_c.B, capacity_init_c.V, capacity_init_c.W, (nx+1, ny+1))

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
solver = MovingDiffusionUnsteadyDiph2(Fluide_1, Fluide_2, bc_b, ic, Δt, Tend, u0, "BE")

# Solve the problem
solve_MovingDiffusionUnsteadyDiph2!(solver, Fluide_1, Fluide_2, u0, Δt, Tend, nt, bc_b, ic, body, body_c, mesh, t, "BE"; method=Base.:\)

# Animation
animate_solution(solver, mesh, body)