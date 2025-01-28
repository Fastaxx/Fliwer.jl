using Fliwer
using IterativeSolvers

### 2D Test Case : Diphasic Steady Diffusion Equation inside a Disk
# Define the mesh
nx, ny = 160, 160
lx, ly = 4., 4.
x0, y0 = 0., 0.
domain = ((x0, lx), (y0, ly))
mesh = CartesianMesh((nx, ny), (lx, ly), (x0, y0))

# Define the body
radius, center = ly/4, (lx/2, ly/2)
circle = Body((x,y,_=0)->sqrt((x-center[1])^2 + (y-center[2])^2) - radius, (x,y,_)->(x,y), domain, false)
circle_c = Body((x,y,_=0)->-(sqrt((x-center[1])^2 + (y-center[2])^2) - radius), (x,y,_)->(x,y), domain, false)

# Identify cells
identify!(mesh, circle)

# Define the capacity
capacity = Capacity(circle, mesh)
capacity_c = Capacity(circle_c, mesh)

# Define the operators
operator = DiffusionOps(capacity.A, capacity.B, capacity.V, capacity.W, (nx+1, ny+1))
operator_c = DiffusionOps(capacity_c.A, capacity_c.B, capacity_c.V, capacity_c.W, (nx+1, ny+1))

# Define the boundary conditions
bc = Dirichlet(1.0)
bc1 = Dirichlet(0.0)
bc2 = Dirichlet(2.0)
bc_b = BorderConditions(Dict{Symbol, AbstractBoundary}(:left => bc1, :right => bc1, :top => bc1, :bottom => bc1))

ic = InterfaceConditions(ScalarJump(1.0, 1.0, 0.0), FluxJump(1.0, 1.0, 0.0))

# Define the source term
f1 = (x,y,_)->4.0 #cos(x)*sin(10*y)
f2 = (x,y,_)->16/10 * (sqrt((x-center[1])^2 + (y-center[2])^2))

# Define the phases
Fluide_1 = Phase(capacity, operator, f1, 1.0)
Fluide_2 = Phase(capacity_c, operator_c, f2, 1.0)

# Define the solver
solver = DiffusionSteadyDiph(Fluide_1, Fluide_2, bc_b, ic)

# Solve the problem
Fliwer.solve_DiffusionSteadyDiph!(solver, Fluide_1, Fluide_2; method=Base.:\)
#solve_DiffusionSteadyDiph!(solver, Fluide_1, Fluide_2; method=IterativeSolvers.gmres, verbose=false)

# Plot the solution usign Makie
plot_solution(solver, mesh, circle, capacity)

# Write the solution to a VTK file
#write_vtk("poisson_2d", mesh, solver)

# Analytical solution
u_analytical(x,y) = 1 - (x-center[1])^2 - (y-center[2])^2

cell_centroids = capacity.C_ω
cell_centroids_c = capacity_c.C_ω
u_ana = map(c -> u_analytical(c[1], c[2]), cell_centroids)
u_ana_c = map(c -> u_analytical(c[1], c[2]), cell_centroids_c)

u_num1 = solver.x[1:end÷2]
u_num2 = solver.x[end÷2+1:end]

u_num1ₒ = u_num1[1:end÷2]
u_num1ᵧ = u_num1[end÷2+1:end]
u_num2ₒ = u_num2[1:end÷2]
u_num2ᵧ = u_num2[end÷2+1:end]

u_ana[capacity.cell_types .== 0] .= NaN
u_ana_c[capacity_c.cell_types .== 0] .= NaN
u_ana = reshape(u_ana, (nx+1, ny+1))
u_ana_c = reshape(u_ana_c, (nx+1, ny+1))

u_num1ₒ[capacity.cell_types .== 0] .= NaN
u_num1ₒ = reshape(u_num1ₒ, (nx+1, ny+1))

u_num2ₒ[capacity_c.cell_types .== 0] .= NaN
u_num2ₒ = reshape(u_num2ₒ, (nx+1, ny+1))

using CairoMakie

fig = Figure()
ax = Axis(fig[1, 1], xlabel="x", ylabel="y", title="Phase 1")
ax1 = Axis(fig[1, 2], xlabel="x", ylabel="y", title="Phase 2")
ax2 = Axis(fig[2, 1], xlabel="x", ylabel="y", title="Analytical solution - Phase 1")
ax3 = Axis(fig[2, 2], xlabel="x", ylabel="y", title="Analytical solution - Phase 2")
hm = heatmap!(ax, u_num1ₒ, colormap=:viridis)
hm1 = heatmap!(ax1, u_num2ₒ, colormap=:viridis)
hm2 = heatmap!(ax2, u_ana, colormap=:viridis)
hm3 = heatmap!(ax3, u_ana_c, colormap=:viridis)
Colorbar(fig[1, 3], hm, label="uₙ(x) - Phase 1")
Colorbar(fig[1, 4], hm1, label="uₙ(x) - Phase 2")
Colorbar(fig[2, 3], hm2, label="uₐ(x) - Phase 1")
Colorbar(fig[2, 4], hm3, label="uₐ(x) - Phase 2")
display(fig)

