using Fliwer
using IterativeSolvers

### 2D Test Case : Advection of a vector field in a cylinder
# Define the mesh
nx, ny = 80, 80
lx, ly = 1., 1.
x0, y0 = 0., 0.
domain = ((x0, lx), (y0, ly))
mesh_p = CartesianMesh((nx, ny), (lx, ly), (x0, y0))

mesh_u = CartesianMesh((nx-1, ny), (lx - lx/nx, ly), (x0 + lx/(2*nx), y0))
mesh_v = CartesianMesh((nx, ny-1), (lx, ly - ly/ny), (x0, y0 + ly/(2*ny)))

@show mesh_p.centers
@show mesh_u.centers
@show mesh_v.centers
# Define the body
radius, center = ly/4, (lx/2, ly/2) .+ (0.01, 0.01)
body = Body((x,y,_=0)->(sqrt((x-center[1])^2 + (y-center[2])^2) - radius), (x,y)->(x,y), domain, true)

# Identify cells in the meshes
identify!(mesh_p, body)
identify!(mesh_u, body)

# Define capacities using the new meshes : Aₐ, Bₐ, Vₐ, Wₐ for a ∈ {p, u, v}
capacity_p = Capacity(body, mesh_p)
capacity_u = Capacity(body, mesh_u)
capacity_v = Capacity(body, mesh_v)

# Initialize the velocity field with a rotating field
uₒx, uₒy = ones(nx*(ny+1)), zeros((nx+1)*ny)
uₒ = (uₒx, uₒy)

# For boundary velocities, if they are zero:
uᵧ = zeros((nx) * (ny+1) + (nx+1) * (ny))

# Initialize the solution vector
x0 = zeros(nx*(ny+1) + (nx+1)*ny + (nx) * (ny+1) + (nx+1) * (ny))

x0[1:nx*(ny+1)] .= uₒx
x0[nx*(ny+1)+1:nx*(ny+1)+(nx+1)*ny] .= uᵧ[1:(nx+1)*ny]
x0[nx*(ny+1)+(nx+1)*ny+1:nx*(ny+1)+(nx+1)*ny+(nx) * (ny+1)] .= uₒy
x0[nx*(ny+1)+(nx+1)*ny+(nx) * (ny+1)+1:end] .= uᵧ[(nx+1)*ny+1:end]

# Define the operators
operator = AdvectionVecOps(capacity_u, capacity_v, (nx+1, ny), uₒ, uᵧ)

# Define the boundary conditions
bc = BorderConditions(Dict{Symbol, AbstractBoundary}(:left => Dirichlet(0.0), :right => Dirichlet(0.0), :top => Dirichlet(0.0), :bottom => Dirichlet(0.0)))

# Define the source term
f = (x, y, z, t) -> 0.0

# Define the phase
velocity = Velocity{2}((capacity_u, capacity_v), operator, f, 1.0, 1.0)

# Define the solver
Δt = 0.1
Tend = 1.0
s = AdvectionVecUnsteadyMono(velocity, bc, Δt, Tend, x0)

# Solve the problem
s.x = IterativeSolvers.gmres(s.A, s.b)

# Plot
using CairoMakie

uₒx = s.x[1:nx*(ny+1)]
uᵧx = s.x[nx*(ny+1)+1:nx*(ny+1)+(nx+1)*ny]
uₒy = s.x[nx*(ny+1)+(nx+1)*ny+1:nx*(ny+1)+(nx+1)*ny+(nx) * (ny+1)]
uᵧy = s.x[nx*(ny+1)+(nx+1)*ny+(nx) * (ny+1)+1:end]

uₒx = reshape(uₒx, nx, ny+1)
uᵧx = reshape(uᵧx, nx+1, ny)
uₒy = reshape(uₒy, nx, ny+1)
uᵧy = reshape(uᵧy, nx+1, ny)

fig = Figure(size = (800, 800))
ax = Axis(fig[1, 1], aspect = DataAspect(), xlabel = "x", ylabel = "y")
heatmap!(ax, uₒx, colormap = :viridis)
Colorbar(fig[1, 2], label = "uₒx")

ax = Axis(fig[2, 1], aspect = DataAspect(), xlabel = "x", ylabel = "y")
heatmap!(ax, uᵧx, colormap = :viridis)
Colorbar(fig[2, 2], label = "uᵧx")

ax = Axis(fig[3, 1], aspect = DataAspect(), xlabel = "x", ylabel = "y")
heatmap!(ax, uₒy, colormap = :viridis)
Colorbar(fig[3, 2], label = "uₒy")

ax = Axis(fig[4, 1], aspect = DataAspect(), xlabel = "x", ylabel = "y")
heatmap!(ax, uᵧy, colormap = :viridis)
Colorbar(fig[4, 2], label = "uᵧy")
display(fig)