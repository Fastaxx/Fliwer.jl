using Fliwer
using IterativeSolvers

### 2D Test Case : Diphasic Unsteady Diffusion Equation with a Disk
# Define the mesh
nx, ny = 40, 40
lx, ly = 4., 4.
x0, y0 = 0., 0.
domain = ((x0, lx), (y0, ly))
mesh = CartesianMesh((nx, ny), (lx, ly), (x0, y0))

# Define the body
radius, center = ly/4, (lx/2, ly/2) .+ (0.01, 0.01)
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
bc = Dirichlet(0.0)
bc_b = BorderConditions(Dict{Symbol, AbstractBoundary}(:left => bc, :right => bc, :top => bc, :bottom => bc))

ic = InterfaceConditions(ScalarJump(1.0, 0.1, 0.0), FluxJump(1.0, 1.0, 0.0))

# Define the source term
f1 = (x,y,z,t)->0.0
f2 = (x,y,z,t)->0.0

# Define the phases
Fluide_1 = Phase(capacity, operator, f1, 10.0)
Fluide_2 = Phase(capacity_c, operator_c, f2, 1.0)

# Initial condition
u0ₒ1 = ones((nx+1)*(ny+1))
u0ᵧ1 = zeros((nx+1)*(ny+1))
u0ₒ2 = zeros((nx+1)*(ny+1))
u0ᵧ2 = zeros((nx+1)*(ny+1))
u0 = vcat(u0ₒ1, u0ᵧ1, u0ₒ2, u0ᵧ2)

# Define the solver
Δt = 0.01
Tend = 1.0
solver = DiffusionUnsteadyDiph(Fluide_1, Fluide_2, bc_b, ic, Δt, Tend, u0)

# Solve the problem
solve!(solver, Fluide_1, Fluide_2, u0, Δt, Tend, bc_b, ic; method=IterativeSolvers.gmres, abstol=1e-15, restart=20, verbose=false)

# Plot the solution using Makie
using CairoMakie

# Reshaper la solution
u1ₒ = reshape(solver.x[1:length(solver.x) ÷ 4], (nx + 1, ny + 1))'
u1ᵧ = reshape(solver.x[length(solver.x) ÷ 4 + 1:2*length(solver.x) ÷ 4], (nx + 1, ny + 1))'
u2ₒ = reshape(solver.x[2*length(solver.x) ÷ 4 + 1:3*length(solver.x) ÷ 4], (nx + 1, ny + 1))'
u2ᵧ = reshape(solver.x[3*length(solver.x) ÷ 4 + 1:end], (nx + 1, ny + 1))'

# Créer les coordonnées x et y
x = range(x0, stop=lx, length=nx+1)
y = range(y0, stop=ly, length=ny+1)

# Tracer la solution en 3D
fig = Figure()

# Phase 1 - Bulk
ax1 = Axis3(fig[1, 1], title = "Phase 1 - Bulk", xlabel = "x", ylabel = "y", zlabel = "u1ₒ")
surface!(ax1, x, y, u1ₒ, colormap = :viridis)

# Phase 1 - Interface
ax2 = Axis3(fig[1, 2], title = "Phase 1 - Interface", xlabel = "x", ylabel = "y", zlabel = "u1ᵧ")
surface!(ax2, x, y, u1ᵧ, colormap = :viridis)

# Phase 2 - Bulk
ax3 = Axis3(fig[2, 1], title = "Phase 2 - Bulk", xlabel = "x", ylabel = "y", zlabel = "u2ₒ")
surface!(ax3, x, y, u2ₒ, colormap = :viridis)

# Phase 2 - Interface
ax4 = Axis3(fig[2, 2], title = "Phase 2 - Interface", xlabel = "x", ylabel = "y", zlabel = "u2ᵧ")
surface!(ax4, x, y, u2ᵧ, colormap = :viridis)

# Afficher le graphique
display(fig)

# Reshaper la solution
u1ₒ = reshape(solver.x[1:length(solver.x) ÷ 4], (nx + 1, ny + 1))'
u1ᵧ = reshape(solver.x[length(solver.x) ÷ 4 + 1:2*length(solver.x) ÷ 4], (nx + 1, ny + 1))'
u2ₒ = reshape(solver.x[2*length(solver.x) ÷ 4 + 1:3*length(solver.x) ÷ 4], (nx + 1, ny + 1))'
u2ᵧ = reshape(solver.x[3*length(solver.x) ÷ 4 + 1:end], (nx + 1, ny + 1))'

# Créer les coordonnées x et y
x = range(x0, stop=lx, length=nx+1)
y = range(y0, stop=ly, length=ny+1)

# Déterminer les limites de la couleur à partir des états
min_val = minimum([minimum(reshape(state[1:length(state) ÷ 4], (nx + 1, ny + 1))') for state in solver.states])
max_val = maximum([maximum(reshape(state[1:length(state) ÷ 4], (nx + 1, ny + 1))') for state in solver.states])

# Créer une figure pour l'animation
fig = Figure(size = (800, 600))

# Phase 1 - Bulk
ax1 = Axis3(fig[1, 1], title = "Phase 1 - Bulk", xlabel = "x", ylabel = "y", zlabel = "u1ₒ")
s1 = surface!(ax1, x, y, u1ₒ, colormap = :viridis, colorrange = (min_val, max_val))

# Phase 1 - Interface
ax2 = Axis3(fig[1, 2], title = "Phase 1 - Interface", xlabel = "x", ylabel = "y", zlabel = "u1ᵧ")
s2 = surface!(ax2, x, y, u1ᵧ, colormap = :viridis, colorrange = (min_val, max_val))

# Phase 2 - Bulk
ax3 = Axis3(fig[2, 1], title = "Phase 2 - Bulk", xlabel = "x", ylabel = "y", zlabel = "u2ₒ")
s3 = surface!(ax3, x, y, u2ₒ, colormap = :viridis, colorrange = (min_val, max_val))

# Phase 2 - Interface
ax4 = Axis3(fig[2, 2], title = "Phase 2 - Interface", xlabel = "x", ylabel = "y", zlabel = "u2ᵧ")
s4 = surface!(ax4, x, y, u2ᵧ, colormap = :viridis, colorrange = (min_val, max_val))

# Fonction pour mettre à jour les surfaces
function update_surfaces!(frame)
    s1[:z] = reshape(solver.states[frame][1:length(solver.states[frame]) ÷ 4], (nx + 1, ny + 1))'
    s2[:z] = reshape(solver.states[frame][length(solver.states[frame]) ÷ 4 + 1:2*length(solver.states[frame]) ÷ 4], (nx + 1, ny + 1))'
    s3[:z] = reshape(solver.states[frame][2*length(solver.states[frame]) ÷ 4 + 1:3*length(solver.states[frame]) ÷ 4], (nx + 1, ny + 1))'
    s4[:z] = reshape(solver.states[frame][3*length(solver.states[frame]) ÷ 4 + 1:end], (nx + 1, ny + 1))'
end

# Créer l'animation
record(fig, "heat_DiffUnsteadyDiph.mp4", 1:length(solver.states); framerate = 10) do frame
    update_surfaces!(frame)
end

# Afficher la figure finale
display(fig)
