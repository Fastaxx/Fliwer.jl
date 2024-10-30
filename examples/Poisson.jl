using Fliwer
using IterativeSolvers

### 2D Test Case : Monophasic Steady Diffusion Equation inside a Disk
# Define the mesh
nx, ny = 80, 80
lx, ly = 4., 4.
x0, y0 = 0., 0.
domain = ((x0, lx), (y0, ly))
mesh = CartesianMesh((nx, ny), (lx, ly), (x0, y0))

# Define the body
radius, center = ly/4, (lx/2, ly/2) .+ (0.01, 0.01)
circle = Body((x,y,_=0)->sqrt((x-center[1])^2 + (y-center[2])^2) - radius, (x,y,_)->(x,y), domain, false)

# Identify cells
identify!(mesh, circle)

# Define the capacity
capacity = Capacity(circle, mesh)

# Define the operators
operator = DiffusionOps(capacity.A, capacity.B, capacity.V, capacity.W, (nx+1, ny+1))

# Define the boundary conditions 
bc = Dirichlet(0.0)

bc_b = BorderConditions(Dict{Symbol, AbstractBoundary}(:left => bc, :right => bc, :top => bc, :bottom => bc))

# Define the source term
f = (x,y,_)->sin(x)*cos(10*y)

Fluide = Phase(capacity, operator, f, 1.0)

# Define the solver
solver = DiffusionSteadyMono(Fluide, bc_b, bc)

# Solve the problem
solve!(solver, Fluide; method=IterativeSolvers.cg, abstol=1e-15, maxiter=1000, verbose=true)

# Plot the solution usign Makie
using CairoMakie
# Reshaper la solution
reshaped_u = reshape(solver.x[1:length(solver.x) ÷ 2], (nx + 1, ny + 1))'

# Tracer la solution avec heatmap
fig = Figure()
ax = Axis(fig[1, 1], title = "Solution Plot", xlabel = "x", ylabel = "y")
hm = heatmap!(ax, mesh.centers[1], mesh.centers[2], reshaped_u, colormap = :viridis)

# Ajouter le zéro de la fonction distance signée de Body
contour!(ax, mesh.nodes[1], mesh.nodes[2], [circle.sdf(xi, yi, 0.0) for yi in mesh.nodes[2], xi in mesh.nodes[1]], levels = [0.0], color = :red, linewidth = 2, label = "Contour SDF=0")

# Ajouter une colorbar
Colorbar(fig[1, 2], hm, label = "Intensity")

# Afficher le graphique
display(fig)

# Fonction pour trouver la valeur maximale d'un tableau
println(maximum(solver.x))