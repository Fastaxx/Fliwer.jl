using Fliwer
using IterativeSolvers
using LinearAlgebra

### 2D Test Case : Monophasic Unsteady Advection-Diffusion Equation inside a Disk
# Define the mesh
nx, ny = 160, 160
lx, ly = 16., 16.
x0, y0 = 0., 0.
domain = ((x0, lx), (y0, ly))
mesh = CartesianMesh((nx, ny), (lx, ly), (x0, y0))

# Define the body : Planar
radius, center = ly/8, (lx/2, ly/2) .+ (0.01, 0.01)
circle = Body((x,y,_=0)->(x - center[2]), (x,y,_)->(x,y), domain, false)

# Identify cells
identify!(mesh, circle)

# Define the capacity
capacity = Capacity(circle, mesh)

# Initialize the velocity field with a Poiseuille flow
# Number of nodes
N = (nx + 1) * (ny + 1)

# Initialize velocity components
uₒx = zeros(N)
uₒy = zeros(N)

# Calculate Poiseuille velocity components
for j in 0:ny
    for i in 0:nx
        idx = i + j * (nx + 1) + 1
        y = (y0 - ly/2) + j * ly/ny
        uₒx[idx] =  0.0
        uₒy[idx] =  y * (1 - y)   # No velocity in the y direction for Poiseuille flow
    end
end

# Interleave the components
uₒ = Vector{Float64}(undef, 2N)
uₒ[1:2:2N] .= uₒx
uₒ[2:2:2N] .= uₒy

# For boundary velocities, if they are zero:
uᵧ = zeros(2N)

# Define the operators
operator = ConvectionOps(capacity.A, capacity.B, capacity.V, capacity.W, (nx+1, ny+1), uₒ, uᵧ)

# Define the boundary conditions
ic = Dirichlet(1.0)
bc = Dirichlet(0.0)

bc_b = BorderConditions(Dict{Symbol, AbstractBoundary}(:left => bc, :right => bc, :top => bc, :bottom => bc))

# Define the source term
f = (x,y,z,t)-> 0.0 #sin(x)*cos(10*y)

Fluide = Phase(capacity, operator, f, 1.0)

# Initial condition
# Analytical solution T (x, t) = 10/(4D (t + 1/2)) exp(− ‖x−xc(t)‖^2/(4D(t+ 1/2)
T0ₒ = zeros((nx+1)*(ny+1))
T0ᵧ = zeros((nx+1)*(ny+1))

# Combine the temperature arrays if needed
T0 = vcat(T0ₒ, T0ᵧ)

# Define the solver
Δt = 0.01
Tend = 2.0
solver = AdvectionDiffusionUnsteadyMono(Fluide, bc_b, ic, Δt, Tend, T0)

# Solve the problem
solve_AdvectionDiffusionUnsteadyMono!(solver, Fluide, T0, Δt, Tend, bc_b, ic; method=IterativeSolvers.bicgstabl, verbose=false, reltol=1e-20)

# Plot the solution
#plot_solution(solver, mesh, circle, capacity)

# Write the solution to a VTK file
write_vtk("advdiff", mesh, solver)

# Animation
#animate_solution(solver, mesh, circle)