using Fliwer
using IterativeSolvers
using SparseArrays

### 1D Test Case : Navier-Stokes
# Define the mesh
nx = 80
lx = 1.
x0 = 0.
domain = ((x0, lx),)
mesh_p = CartesianMesh((nx,), (lx,), (x0,))

mesh_u = CartesianMesh((nx-1,), (lx - lx/nx,), (x0 + lx/(2*nx),))

# Define the body
body = NoBody1D(domain)

# Identify cells in the meshes
identify!(mesh_p, body)
identify!(mesh_u, body)

# Define capacities using the new meshes : Aₐ, Bₐ, Vₐ, Wₐ for a ∈ {p, u}
capacity_p = Capacity(body, mesh_p)
capacity_u = Capacity(body, mesh_u)

# Initialize the velocity field with a divergence-free field : 
uₒ = (ones(nx),)
uᵧ = zeros(nx)

# Initialize the solution vector
x0 = zeros(4nx)
x0[1:nx] .= uₒ[1] # Bulk velocity uₒ
x0[nx+1:2nx] .= uᵧ # Boundary velocity uᵧ
x0[2nx+1:3nx] .= zeros(nx) # Pressure pₒ
x0[3nx+1:4nx] .= zeros(nx) # Pressure pᵧ

# Define the operators
operator = NavierStokesOps(capacity_u, (nx,), uₒ, uᵧ)

# Define the boundary conditions
bc = BorderConditions(Dict{Symbol, AbstractBoundary}(:top => Dirichlet(0.0), :bottom => Dirichlet(0.0)))

# Define the source term
f = (x, y, z, t) -> 0.0

# Define the phase
ρ = 1.0
Re = 1.0
#pressure = Pressure(capacity_p, operator, f, 1.0)
velocity = Velocity{1}((capacity_u,), operator, f, ρ, Re)

# Define the solver
Δt = 0.001
Tend = 0.0
solver = NavierStokesUnsteadyMono(velocity, bc, Δt, Tend, x0)

# Solve the problem
solve_NavierStokesUnsteadyMono!(solver, velocity, Δt, Tend, bc; method=IterativeSolvers.bicgstabl, abstol=1e-15, verbose=false)

# Plot
using CairoMakie

uₒ = solver.x[1:nx]
uᵧ = solver.x[nx+1:2nx]
pₒ = solver.x[2nx+1:3nx]
pᵧ = solver.x[3nx+1:4nx]

fig = Figure(size = (800, 400))
ax1 = Axis(fig, xlabel = "x", ylabel = "Pressure", title = "Pressure field")
ax2 = Axis(fig, xlabel = "x", ylabel = "Velocity", title = "Velocity field")
ax3 = Axis(fig, xlabel = "x", ylabel = "Velocity", title = "Velocity field")
ax4 = Axis(fig, xlabel = "x", ylabel = "Pressure", title = "Pressure field")

lines!(ax1, pₒ, color = :blue, linewidth = 2)
lines!(ax2, uₒ, color = :red, linewidth = 2, label = "uₒ")
lines!(ax3, uᵧ, color = :green, linewidth = 2, label = "uᵧ")
lines!(ax4, pᵧ, color = :blue, linewidth = 2)


fig[1, 1] = ax1
fig[1, 2] = ax2
fig[1, 3] = ax3
fig[1, 4] = ax4

display(fig)
