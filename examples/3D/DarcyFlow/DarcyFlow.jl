using Fliwer
using IterativeSolvers

### 3D Test Case : Darcy Flow with a Sphere
# Define the mesh
nx, ny, nz = 40, 40, 40
lx, ly, lz = 20., 10., 10.
x0, y0, z0 = -20., -10., -10.
domain = ((x0, lx), (y0, ly), (z0, lz))
mesh = CartesianMesh((nx, ny, nz), (lx, ly, lz), (x0, y0, z0))

# three new staggered grids are defined, each corresponding to a velocity component, one displaced half a cell width, the other one displaced half a cell height and the last one displaced half a cell depth.
mesh_u = CartesianMesh((nx-1, ny, nz), (lx - lx/nx, ly, lz), (x0 + lx/(2*nx), y0, z0))
mesh_v = CartesianMesh((nx, ny-1, nz), (lx, ly - ly/ny, lz), (x0, y0 + ly/(2*ny), z0))
mesh_w = CartesianMesh((nx, ny, nz-1), (lx, ly, lz - lz/nz), (x0, y0, z0 + lz/(2*nz)))

# Define the body
radius, center = ly/2, (lx/2, ly/2, lz/2) #.+ (0.01, 0.01)
sphere = Body((x,y,z)->-(sqrt((x-center[1])^2 + (y-center[2])^2 + (z-center[3])^2) - radius), (x,y,z,_)->(x,y,z), domain, false)

# Identify cells
identify!(mesh, sphere)
identify!(mesh_u, sphere)
identify!(mesh_v, sphere)
identify!(mesh_w, sphere)

# Define the capacity
capacity = Capacity(sphere, mesh)
capacity_u = Capacity(sphere, mesh_u)
capacity_v = Capacity(sphere, mesh_v)
capacity_w = Capacity(sphere, mesh_w)

# Define the operators
operator_p = DiffusionOps(capacity.A, capacity.B, capacity.V, capacity.W, (nx+1, ny+1, nz+1))
operator_u = DiffusionOps(capacity_u.A, capacity_u.B, capacity_u.V, capacity_u.W, (nx, ny+1, nz+1))
operator_v = DiffusionOps(capacity_v.A, capacity_v.B, capacity_v.V, capacity_v.W, (nx+1, ny, nz+1))
operator_w = DiffusionOps(capacity_w.A, capacity_w.B, capacity_w.V, capacity_w.W, (nx+1, ny+1, nz))

# Define the boundary conditions for pressure for left and right faces. Others faces don't have BC
bc_10 = Dirichlet(10.0)
bc_20 = Dirichlet(20.0)

bc_p = BorderConditions(Dict{Symbol, AbstractBoundary}(:top => bc_10, :bottom => bc_20))

ic = Dirichlet(0.0)

# Define the source term
f = (x,y,z)-> 0.0

# Define the phase
K = 1.0
Fluide = Phase(capacity, operator_p, f, K)

# Define the solver
solver = DarcyFlow(Fluide, bc_p, ic)

# Solve the problem
solve_DarcyFlow!(solver, Fluide; method=Base.:\)

# Write the solution to a VTK file
write_vtk("darcy_3d", mesh, solver)

plot_solution(solver, mesh, sphere, capacity)

function ∇(operator::AbstractOperators, p::Vector{Float64})
    ∇ = operator.Wꜝ * (operator.G * p[1:div(end,2)] + operator.H * p[div(end,2)+1:end])

    return ∇
end

# Compute the velocity field from the pressure field
u = ∇(operator_p, solver.x)

ux = u[1:div(end,3)]
uy = u[div(end,3)+1:2*div(end,3)]
uz = u[2*div(end,3)+1:end]

ux = reshape(ux, (nx+1, ny+1, nz+1))

"""

# Plot the velocity field
using CairoMakie

fig = Figure()
ax = Axis3(fig,  xlabel = "x", ylabel = "y", zlabel = "z", aspect = :equal)

heatmap!(ax, ux[:,:,div(end,3)], colormap = :viridis)

fig
"""