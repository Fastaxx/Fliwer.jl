using Fliwer
using IterativeSolvers

### 1D Test Case : Monophasic Unsteady Diffusion Equation inside a moving body
# Define the spatial mesh
nx = 10
lx = 4.
x0 = 0.
domain = ((x0, lx),)
mesh = CartesianMesh((nx,), (lx,), (x0,))

# Define the time mesh
Δt = 0.01
Tend = 1.0
nt = Int(Tend/Δt)
t = [i*Δt for i in 0:nt]
t = t[2:3]

# Define the body
xf = lx/4   # Interface position
c = 0.1     # Interface velocity
body = Body((x,t, _=0)->(x - xf - c*t), (x,)->(x,), domain, false)  # Body moving to the right

# Define the space-time mesh
spaceTimeMesh = CartesianSpaceTimeMesh(mesh, t)

# Identify cells TO FIX
#identify!(spaceTimeMesh, body)

# Define the capacity
capacity = Capacity(body, spaceTimeMesh)

# Define the operators
operator = SpaceTimeOps(capacity.A, capacity.B, capacity.V, capacity.W, (nx+1, 2))

# Define the boundary conditions
bc = Dirichlet(0.0)
bc1 = Dirichlet(0.0)

bc_b = BorderConditions(Dict{Symbol, AbstractBoundary}(:left => bc, :right => bc))

# Define the source term
f = (x,y,z,t)-> 4.0 #sin(x)*cos(10*y)

Fluide = Phase(capacity, operator, f, 1.0)

# Initial condition
u0ₒ = zeros((nx+1))
u0ᵧ = zeros((nx+1))
u0 = vcat(u0ₒ, u0ᵧ)

# Define the solver
