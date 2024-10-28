using Fliwer

### 2D Test Case : Diphasic Steady Diffusion Equation inside a Disk
# Define the mesh
nx, ny = 80, 80
lx, ly = 4., 4.
x0, y0 = 0., 0.
domain = ((x0, lx), (y0, ly))
mesh = CartesianMesh((nx, ny), (lx, ly), (x0, y0))

# Define the body
radius, center = ly/4, (lx/2, ly/2) .+ (0.01, 0.01)
circle = Body((x,y,_=0)->sqrt((x-center[1])^2 + (y-center[2])^2) - radius, (x,y,_)->(x,y), domain, false)
circle_c = Body((x,y,_=0)->-(sqrt((x-center[1])^2 + (y-center[2])^2) - radius), (x,y,_)->(x,y), domain, false)

# Define the capacity
capacity = Capacity(circle, mesh)
capacity_c = Capacity(circle_c, mesh)

# Define the operators
operator = DiffusionOps(capacity.A, capacity.B, capacity.V, capacity.W, (nx+1, ny+1))
operator_c = DiffusionOps(capacity_c.A, capacity_c.B, capacity_c.V, capacity_c.W, (nx+1, ny+1))

# Define the boundary conditions
bc = Dirichlet(0.0)
bc_b = BorderConditions(Dict{Symbol, AbstractBoundary}(:left => bc, :right => bc, :top => bc, :bottom => bc))

ic = InterfaceConditions(ScalarJump(1.0, 0.0, 0.0), FluxJump(1.0, 1.0, 0.0))

# Define the source term
f1 = (x,y,_)->1.0
f2 = (x,y,_)->0.0

# Define the phases
Fluide_1 = Phase(capacity, operator, f1, 1.0)
Fluide_2 = Phase(capacity_c, operator_c, f2, 1.0)

# Define the solver
solver = DiffusionSteadyDiph(Fluide_1, Fluide_2, bc_b, ic)

# Solve the problem
u = solve!(solver, Fluide_1, Fluide_2)

# Plot the solution usign Makie

using CairoMakie

# Reshaper la solution
u1ₒ = reshape(u[1:length(u) ÷ 4], (nx + 1, ny + 1))'
u1ᵧ = reshape(u[length(u) ÷ 4 + 1:2*length(u) ÷ 4], (nx + 1, ny + 1))'
u2ₒ = reshape(u[2*length(u) ÷ 4 + 1:3*length(u) ÷ 4], (nx + 1, ny + 1))'
u2ᵧ = reshape(u[3*length(u) ÷ 4 + 1:end], (nx + 1, ny + 1))'

# Tracer la solution avec heatmap
fig = Figure()

# Phase 1 - Bulk
ax1 = Axis(fig[1, 1], title = "Phase 1 - Bulk", xlabel = "x", ylabel = "y")
hm1 = heatmap!(ax1, u1ₒ, colormap = :viridis)
cb1 = Colorbar(fig[1, 2], hm1, label = "Intensity")

# Phase 1 - Interface
ax2 = Axis(fig[1, 3], title = "Phase 1 - Interface", xlabel = "x", ylabel = "y")
hm2 = heatmap!(ax2, u1ᵧ, colormap = :viridis)
cb2 = Colorbar(fig[1, 4], hm2, label = "Intensity")

# Phase 2 - Bulk
ax3 = Axis(fig[2, 1], title = "Phase 2 - Bulk", xlabel = "x", ylabel = "y")
hm3 = heatmap!(ax3, u2ₒ, colormap = :viridis)
cb3 = Colorbar(fig[2, 2], hm3, label = "Intensity")

# Phase 2 - Interface
ax4 = Axis(fig[2, 3], title = "Phase 2 - Interface", xlabel = "x", ylabel = "y")
hm4 = heatmap!(ax4, u2ᵧ, colormap = :viridis)
cb4 = Colorbar(fig[2, 4], hm4, label = "Intensity")

# Afficher le graphique
display(fig)
