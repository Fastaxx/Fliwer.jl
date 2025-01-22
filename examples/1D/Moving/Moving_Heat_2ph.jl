using Fliwer
using IterativeSolvers
using LinearAlgebra
using SparseArrays

### 1D Test Case : Diphasic Unsteady Diffusion Equation inside a moving body
# Define the spatial mesh
nx = 40
lx = 10.
x0 = 0.
domain = ((x0, lx),)
mesh = CartesianMesh((nx,), (lx,), (x0,))

# Define the time mesh
Δt = 0.001
Tend = 0.1
nt = Int(Tend/Δt)
t = [i*Δt for i in 0:nt]

# Define the body
xf = 0.53*lx   # Interface position
c = 1.0     # Interface velocity
initial_body = Body((x,_=0)->(x - xf), (x,_)->(x), domain, false)  # Initial body
initial_body_c = Body((x,_=0)->-(x - xf), (x,_)->(x), domain, false)  # Initial body complement

body = Body((x,t, _=0)->(x - xf - c*sqrt(t)), (x,)->(x,), domain, false)  # Body moving to the right
body_c = Body((x,t, _=0)->-(x - xf - c*sqrt(t)), (x,)->(x,), domain, false)  # Body moving to the left

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
operator = SpaceTimeOps(capacity.A, capacity.B, capacity.V, capacity.W, (nx+1, 2))
operator_c = SpaceTimeOps(capacity_c.A, capacity_c.B, capacity_c.V, capacity_c.W, (nx+1, 2))

operator_init = DiffusionOps(capacity_init.A, capacity_init.B, capacity_init.V, capacity_init.W, (nx+1,))
operator_init_c = DiffusionOps(capacity_init_c.A, capacity_init_c.B, capacity_init_c.V, capacity_init_c.W, (nx+1,))

# Define the boundary conditions
bc = Dirichlet(1.0)
bc1 = Dirichlet(1.0)
bc_b = BorderConditions(Dict{Symbol, AbstractBoundary}(:top => bc, :bottom => bc1))

ic = InterfaceConditions(ScalarJump(1.0, 1.0, 0.0), FluxJump(1.0, 1.0, 0.0))

# Define the source term
f1 = (x,y,z,t)->0.0
f2 = (x,y,z,t)->0.0

Fluide_1 = Phase(capacity, operator, f1, 1.0)
Fluide_2 = Phase(capacity_c, operator_c, f2, 1.0)

# Initial condition
u0ₒ1 = ones(nx+1)
u0ᵧ1 = ones(nx+1)
u0ₒ2 = ones(nx+1)
u0ᵧ2 = ones(nx+1)

u0 = vcat(u0ₒ1, u0ᵧ1, u0ₒ2, u0ᵧ2)

# Define the solver
solver = MovingDiffusionUnsteadyDiph(Fluide_1, Fluide_2, bc_b, ic, Δt, Tend, u0, "BE")

# Solve the problem
solve_MovingDiffusionUnsteadyDiph!(solver, Fluide_1, Fluide_2, u0, Δt, Tend, nt, bc_b, ic, body, mesh, t, "BE"; method=Base.:\)

# Plot the solution
#plot_solution(solver, mesh, body, capacity)

# Animation
animate_solution(solver, mesh, body)