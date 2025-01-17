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
Δt = 0.001
Tend = 0.5
nt = Int(Tend/Δt)
t = [i*Δt for i in 0:nt]

# Define the body
xf = 0.53*lx   # Interface position
c = -1.0     # Interface velocity
initial_body = Body((x,_=0)->(x - xf), (x,_)->(x), domain, false)  # Initial body
body = Body((x,t, _=0)->(x - xf + c*sqrt(t)), (x,)->(x,), domain, false)  # Body moving to the right

function plot_body(body, x0, lx, nx, Tend, nt)
    x = range(x0, stop=lx, length=nx)
    t = range(0, stop=Tend, length=nt)
    fig = Figure()
    ax = Axis(fig[1, 1])
    for i in 1:length(t)
        body_pos = [body.sdf(x, t[i]) for x in x]
        lines!(ax, x, body_pos, color=:blue, linewidth=2, label="t=$i")
    end
    display(fig)
end

# Define the space-time mesh
spaceTimeMesh = CartesianSpaceTimeMesh(mesh, t[2:3])

# Identify cells
identify!(mesh, initial_body)
spaceTimeMesh.tag = mesh.tag

# Define the capacity
capacity = Capacity(body, spaceTimeMesh)
capacity_init = Capacity(initial_body, mesh)

# Define the operators
operator = SpaceTimeOps(capacity.A, capacity.B, capacity.V, capacity.W, (nx+1, 2))
operator_init = DiffusionOps(capacity_init.A, capacity_init.B, capacity_init.V, capacity_init.W, (nx+1,))


A = operator.Wꜝ*(operator.G*ones(size(operator.G,2)) + operator.H*ones(size(operator.H,2)))
@show A

Aω = A[1:end÷2]
Aᵧ = A[end÷2+1:end]

Aω = reshape(Aω, (nx+1, 2))
Aᵧ = reshape(Aᵧ, (nx+1, 2))

using CairoMakie

fig = Figure()
ax = Axis(fig[1, 1])
hm = heatmap!(ax, Aω, colormap=:viridis)
Colorbar(fig[1, 2], hm, label="Aω")
display(fig)

readline()

# Define the boundary conditions
bc = Dirichlet(1.0)
bc1 = Dirichlet(1.0)

bc_b = BorderConditions(Dict{Symbol, AbstractBoundary}(:top => Dirichlet(0.0), :bottom => bc1))

# Define the source term
f = (x,y,z,t)-> 0.0 #sin(x)*cos(10*y)

Fluide = Phase(capacity, operator, f, 1.0)

# Initial condition
u0ₒ = ones((nx+1))
u0ᵧ = ones((nx+1))
u0 = vcat(u0ₒ, u0ᵧ)

# Define the solver
solver = MovingDiffusionUnsteadyMono(Fluide, bc_b, bc, Δt, Tend, u0)

# Solve the problem
solve_MovingDiffusionUnsteadyMono!(solver, Fluide, u0, Δt, Tend, nt, bc_b, bc, body, mesh, t; method=Base.:\)

# Write the solution to a VTK file
#write_vtk("moving_heat_1d", mesh, solver)

# Plot the solution
plot_solution(solver, mesh, body, capacity; state_i=1)

# Animation
animate_solution(solver, mesh, body)
