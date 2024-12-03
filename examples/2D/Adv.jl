using Fliwer
using IterativeSolvers

### 2D Test Case : Monophasic Steady Advection Equation inside a Disk
# Define the mesh
nx, ny = 40, 40
lx, ly = 4., 4.
x0, y0 = 0., 0.
domain = ((x0, lx), (y0, ly))
mesh = CartesianMesh((nx, ny), (lx, ly), (x0, y0))

# Define the body
radius, center = ly/4, (lx/2, ly/2) .+ (0.01, 0.01)
circle = Body((x,y,_=0)->(sqrt((x-center[1])^2 + (y-center[2])^2) - radius), (x,y,_)->(x,y), domain, false)

# Identify cells
identify!(mesh, circle)

# Define the capacity
capacity = Capacity(circle, mesh)

# Initialize the velocity field by setting the velocity to zero
uₒx, uₒy = ones((nx+1)*(ny+1)), ones((nx+1)*(ny+1))
uγx, uγy = zeros((nx+1)*(ny+1)), zeros((nx+1)*(ny+1))

uₒ, uᵧ = vcat(uₒx, uₒy), vcat(uγx, uγy)

# Define the operators
operator = ConvectionOps(capacity.A, capacity.B, capacity.V, capacity.W, (nx+1, ny+1), uₒ, uᵧ)

# Define the boundary conditions
ic = Dirichlet(0.0)
bc = Dirichlet(0.0)

bc_b = BorderConditions(Dict{Symbol, AbstractBoundary}(:left => bc, :right => bc, :top => bc, :bottom => bc))

# Define the source term
f = (x,y,z,t)-> 0.0 #sin(x)*cos(10*y)

Fluide = Phase(capacity, operator, f, 1.0)

# Initial condition
u0ₒ = zeros((nx+1)*(ny+1))
u0ᵧ = zeros((nx+1)*(ny+1))
x, y = mesh.centers[1], mesh.centers[2]

for j in 1:ny
    for i in 1:nx
        u0ₒ[i + (j-1)*(nx+1)] = 10/(4*1/2) * exp(-((x[i]-center[1])^2 + (y[j]-center[2])^2)/(4*1/2))
        u0ᵧ[i + (j-1)*(nx+1)] = 10/(4*1/2) * exp(-((x[i]-center[1])^2 + (y[j]-center[2])^2)/(4*1/2))
    end
end

u0 = vcat(u0ₒ, u0ᵧ)

# Define the solver
Δt = 0.01
Tend = 1.0
solver = AdvectionUnsteadyMono(Fluide, bc_b, ic, Δt, Tend, u0)

# Solve the problem
solve_AdvectionUnsteadyMono!(solver, Fluide, u0, Δt, Tend, bc_b, ic; method=IterativeSolvers.gmres, abstol=1e-15, verbose=false)

# Plot the solution
plot_solution(solver, mesh, circle, capacity)

# Write the solution to a VTK file
write_vtk("advection", mesh, solver)