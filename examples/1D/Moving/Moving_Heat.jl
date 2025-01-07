using Fliwer
using IterativeSolvers

### 1D Test Case : Monophasic Unsteady Diffusion Equation inside a moving body
# Define the spatial mesh
nx = 50
lx = 10.
x0 = 0.
domain = ((x0, lx),)
mesh = CartesianMesh((nx,), (lx,), (x0,))

# Define the time mesh
Δt = 0.01
Tend = 1.0
nt = Int(Tend/Δt)
@show nt
t = [i*Δt for i in 0:nt]

# Define the body
xf = lx/4   # Interface position
c = 1.0     # Interface velocity
body = Body((x,t, _=0)->(x - xf - c*t), (x,)->(x,), domain, false)  # Body moving to the right

# Define the space-time mesh
spaceTimeMesh = CartesianSpaceTimeMesh(mesh, t[2:3])

# Identify cells TO FIX
#identify!(spaceTimeMesh, body)

# Define the capacity
capacity = Capacity(body, spaceTimeMesh)

# Define the operators
operator = SpaceTimeOps(capacity.A, capacity.B, capacity.V, capacity.W, (nx+1, 2))

# Define the boundary conditions
bc = Dirichlet(1.0)
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
solver = MovingDiffusionUnsteadyMono(Fluide, bc_b, bc, Δt, Tend, u0)

# Solve the problem
solve_MovingDiffusionUnsteadyMono!(solver, Fluide, u0, Δt, Tend, nt, bc_b, bc, body, mesh, t; method=IterativeSolvers.gmres, restart=10, maxiter=1000, verbose=false)

# Write the solution to a VTK file
write_vtk("moving_heat_1d", mesh, solver)

# Plot the solution
#plot_solution(solver, mesh, body, capacity; state_i=1)

