using Fliwer
using IterativeSolvers
using LinearAlgebra
using SparseArrays

### 1D Test Case : Monophasic Unsteady Diffusion Equation inside a moving body
# Define the spatial mesh
nx = 40
lx = 10.
x0 = 0.
domain = ((x0, lx),)
mesh = CartesianMesh((nx,), (lx,), (x0,))

# Define the time mesh
# Example of adaptive time stepping
function create_adaptive_times(Tend::Float64, nsteps::Int, ratio::Float64=1.1)
    times = zeros(Float64, nsteps+1)
    times[1] = 0.0
    dt = Tend * (1 - ratio) / (1 - ratio^nsteps)
    for i in 2:nsteps+1
        times[i] = times[i-1] + dt * ratio^(i-2)
    end
    return times
end

Δt = 0.001
Tend = 0.2
nt = Int(Tend/Δt)
t = [i*Δt for i in 0:nt]

# Define the body
xf = 0.52*lx   # Interface position
c = 1.0     # Interface velocity
initial_body = Body((x,_=0)->(x - xf), (x,_)->(x), domain, false)  # Initial body
body = Body((x,t, _=0)->(x - xf - c*sqrt(t)), (x,)->(x,), domain, false)  # Body moving to the right

# Define the space-time mesh
spaceTimeMesh = CartesianSpaceTimeMesh(mesh, t[1:2])

# Identify cells
identify!(mesh, initial_body)
spaceTimeMesh.tag = mesh.tag

# Define the capacity
capacity = Capacity(body, spaceTimeMesh)
capacity_init = Capacity(initial_body, mesh)

# Define the operators
operator = SpaceTimeOps(capacity.A, capacity.B, capacity.V, capacity.W, (nx+1, 2))
operator_init = DiffusionOps(capacity_init.A, capacity_init.B, capacity_init.V, capacity_init.W, (nx+1,))

# Define the boundary conditions
bc = Dirichlet(1.0)
bc1 = Dirichlet(0.0)

bc_b = BorderConditions(Dict{Symbol, AbstractBoundary}(:top => Dirichlet(0.0), :bottom => bc1))

# Define the source term
f = (x,y,z,t)-> 0.0 #sin(x)*cos(10*y)

Fluide = Phase(capacity, operator, f, 1.0)

# Initial condition
u0ₒ = zeros((nx+1))
u0ᵧ = ones((nx+1))
u0 = vcat(u0ₒ, u0ᵧ)

# Define the solver
solver = MovingDiffusionUnsteadyMono(Fluide, bc_b, bc, Δt, Tend, u0, "CN")

# Solve the problem
solve_MovingDiffusionUnsteadyMono!(solver, Fluide, u0, Δt, Tend, nt, bc_b, bc, body, mesh, t, "CN"; method=Base.:\)

# Write the solution to a VTK file
#write_vtk("moving_heat_1d", mesh, solver)

# Plot the solution
plot_solution(solver, mesh, body, capacity; state_i=1)

# Animation
animate_solution(solver, mesh, body)
