using Fliwer
using IterativeSolvers

### 1D Test Case : Monophasic Unsteady Advection Equation
# Define the mesh
nx = 160
lx = 4.0
x0 = 0.0
domain=((x0,lx),)
mesh = CartesianMesh((nx,), (lx,), (x0,))

# Define the body
xint = 2.0 + 0.1
body = Body((x, _=0) -> (x - xint),(x,_=0)->(x),domain,false)

# Identify cells
identify!(mesh, body)

# Define the capacity
capacity = Capacity(body, mesh)

# Initialize the velocity field by setting the velocity to zero
uₒ = zeros(nx+1)
uᵧ = zeros(nx+1)

# Define the operators
operator = ConvectionOps(capacity.A, capacity.B, capacity.V, capacity.W, (nx+1,), uₒ, uᵧ)

# Define the boundary conditions
bc0 = Dirichlet(0.0)
bc1 = Dirichlet(0.0)

bc_b = BorderConditions(Dict{Symbol, AbstractBoundary}(:top => bc0, :bottom => bc1))

# Define the source term
f = (x,y,z,t)->0.0

# Define the phase
Fluide = Phase(capacity, operator, f, 1.0)

# Initial condition
u0ₒ = zeros(nx+1)
u0ᵧ = zeros(nx+1)
u0 = vcat(u0ₒ, u0ᵧ)

# Define the solver
Δt = 0.01
Tend = 1.0
solver = AdvectionUnsteadyMono(Fluide, bc_b, bc0, Δt, Tend, u0)

# Solve the problem
solve!(solver, Fluide, u0, Δt, Tend, bc_b, bc0; method=IterativeSolvers.cg, abstol=1e-15, verbose=false)

# Write the solution to a VTK file
write_vtk("advection", mesh, solver)

# Plot the solution
plot_solution(solver, mesh, body, capacity; state_i=10)