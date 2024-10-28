using Fliwer

### 2D Test Case : Monophasic Unsteady Diffusion Equation inside a Disk
# Define the mesh
nx, ny = 40, 40
lx, ly = 4., 4.
x0, y0 = 0., 0.
domain = ((x0, lx), (y0, ly))
mesh = CartesianMesh((nx, ny), (lx, ly), (x0, y0))

# Define the body
radius, center = ly/4, (lx/2, ly/2)
circle = Body((x,y,_=0)->sqrt((x-center[1])^2 + (y-center[2])^2) - radius, (x,y,_)->(x,y), domain, false)

# Define the capacity
capacity = Capacity(circle, mesh)

# Define the operators
operator = DiffusionOps(capacity.A, capacity.B, capacity.V, capacity.W, (nx+1, ny+1))

# Define the boundary conditions 
bc = Dirichlet(10.0)

bc_b = BorderConditions(Dict{Symbol, AbstractBoundary}(:left => bc, :right => bc, :top => bc, :bottom => bc))

# Define the source term
f = (x,y,z,t)->0.0

# Define the phase
Fluide = Phase(capacity, operator, f, 1.0)

# Initial condition
u0ₒ = zeros((nx+1)*(ny+1))
u0ᵧ = ones((nx+1)*(ny+1))
u0 = vcat(u0ₒ, u0ᵧ)

# Define the solver
Δt = 0.001
Tend = 1.0
solver = DiffusionUnsteadyMono(Fluide, bc_b, bc, Δt, Tend, u0)

# Solve the problem
u, states = solve!(solver, Fluide, u0, Δt, Tend, bc_b, bc)

# Plot the solution using Makie
using CairoMakie

# Reshaper la solution
uₒ = reshape(u[1:length(u) ÷ 2], (nx + 1, ny + 1))'
uᵧ = reshape(u[length(u) ÷ 2 + 1:end], (nx + 1, ny + 1))'

# Tracer la solution avec heatmap
fig = Figure()

# Phase 1 - Bulk
ax1 = Axis(fig[1, 1], title = "Phase 1 - Bulk", xlabel = "x", ylabel = "y", aspect = DataAspect())
hm1 = heatmap!(ax1, uₒ, colormap = :viridis)
cb1 = Colorbar(fig[1, 2], hm1, label = "Intensity")

# Phase 2 - Interface
ax2 = Axis(fig[1, 3], title = "Phase 2 - Interface", xlabel = "x", ylabel = "y", aspect = DataAspect())
hm2 = heatmap!(ax2, uᵧ, colormap = :viridis)
cb2 = Colorbar(fig[1, 4], hm2, label = "Intensity")

# Afficher le graphique
display(fig)

# Créer une figure pour l'animation
fig = Figure(size = (800, 400))
ax = Axis(fig[1, 1], title = "Diffusion Unsteady Monophasic", xlabel = "x", ylabel = "y", aspect = DataAspect())

# Fixer les limites des axes
xlims!(ax, 0, nx)
ylims!(ax, 0, ny)

# Déterminer les limites de la couleur à partir des états
min_val = minimum([minimum(reshape(state[1:length(state) ÷ 2], (nx + 1, ny + 1))') for state in states])
max_val = maximum([maximum(reshape(state[1:length(state) ÷ 2], (nx + 1, ny + 1))') for state in states])

# Tracer la première heatmap avec les limites de couleur fixées
hm = heatmap!(ax, reshape(states[1][1:length(states[1]) ÷ 2], (nx + 1, ny + 1))', colormap = :viridis, colorrange = (min_val, max_val))
cb = Colorbar(fig[1, 2], hm, label = "Intensity")

# Fonction pour mettre à jour l'animation
function update_heatmap!(frame)
    hm[1] = reshape(states[frame][1:length(states[frame]) ÷ 2], (nx + 1, ny + 1))'
end

# Créer l'animation
record(fig, "heat_DiffUnsteadyMono.mp4", 1:length(states); framerate = 10) do frame
    update_heatmap!(frame)
end

# Afficher la figure
display(fig)
