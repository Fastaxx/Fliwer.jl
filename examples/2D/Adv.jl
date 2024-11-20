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
circle = Body((x,y,_=0)->-(sqrt((x-center[1])^2 + (y-center[2])^2) - radius), (x,y,_)->(x,y), domain, false)

# Identify cells
identify!(mesh, circle)

# Define the capacity
capacity = Capacity(circle, mesh)

# Initialize the velocity field by setting the velocity to zero
uₒx, uₒy = zeros((nx+1)*(ny+1)), zeros((nx+1)*(ny+1))
uγx, uγy = zeros((nx+1)*(ny+1)), zeros((nx+1)*(ny+1))

uₒ, uᵧ = vcat(uₒx, uₒy), vcat(uγx, uγy)

# Define the operators
operator = ConvectionOps(capacity.A, capacity.B, capacity.V, capacity.W, (nx+1, ny+1), uₒ, uᵧ)

# Define the boundary conditions
bc = Dirichlet(1.0)
bc1 = Dirichlet(0.0)

bc_b = BorderConditions(Dict{Symbol, AbstractBoundary}(:left => bc, :right => bc, :top => bc1, :bottom => bc))

# Define the source term
f = (x,y,_)-> 4.0 #sin(x)*cos(10*y)

Fluide = Phase(capacity, operator, f, 1.0)

# Define the solver
solver = AdvectionSteadyMono(Fluide, bc_b, bc)

# Solve the problem
solve!(solver, Fluide; method=IterativeSolvers.bicgstabl, verbose=false)