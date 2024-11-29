using Fliwer
using IterativeSolvers

### 1D Test Case : Monophasic Unsteady Diffusion Equation 
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

# Define the operators
operator = DiffusionOps(capacity.A, capacity.B, capacity.V, capacity.W, (nx+1,))

# Define the boundary conditions
bc0 = Dirichlet(0.0)
bc1 = Dirichlet(1.0)

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
solver = DiffusionUnsteadyMono(Fluide, bc_b, bc0, Δt, Tend, u0)

# Solve the problem
solve!(solver, Fluide, u0, Δt, Tend, bc_b, bc0; method=IterativeSolvers.gmres, restart=10, maxiter=1000, verbose=false)

# Write the solution to a VTK file
write_vtk("heat_1d", mesh, solver)

# Plot the solution
plot_solution(solver, mesh, body, capacity; state_i=10)

# Animation
animate_solution(solver, mesh, body)

