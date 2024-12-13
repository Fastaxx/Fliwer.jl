using Fliwer
using IterativeSolvers
using SparseArrays
using LinearAlgebra
using CairoMakie

### 2D Test Case : Monophasic Steady Diffusion Equation inside a Disk
# Define the mesh
nx, ny = 80, 80
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

# Define the operators
operator = DiffusionOps(capacity.A, capacity.B, capacity.V, capacity.W, (nx+1, ny+1))

# Define the boundary conditions 
bc = Dirichlet(1.0)
bc1 = Dirichlet(0.0)

bc_b = BorderConditions(Dict{Symbol, AbstractBoundary}(:left => bc1, :right => bc1, :top => bc1, :bottom => bc1))

# Define the source term
f = (x,y,_)-> 4.0 #sin(x)*cos(10*y)

Fluide = Phase(capacity, operator, f, 1.0)

# Define the solver
solver = DiffusionSteadyMono(Fluide, bc_b, bc)

function remove_zero_rows_cols!(A::SparseMatrixCSC{Float64, Int}, b::Vector{Float64})
    # Compute sums of absolute values along rows and columns
    row_sums = vec(sum(abs.(A), dims=2))
    col_sums = vec(sum(abs.(A), dims=1))

    # Find indices of non-zero rows and columns
    rows_idx = findall(row_sums .!= 0.0)
    cols_idx = findall(col_sums .!= 0.0)

    # Create new matrix and RHS vector
    A = A[rows_idx, cols_idx]
    b = b[rows_idx]

    return A, b, rows_idx, cols_idx
end

@time Fliwer.solve!(solver, Fluide; method=IterativeSolvers.bicgstabl, verbose=false, reltol=1e-20)

# Store the size of the original system
n = Int(size(solver.A, 1) /2)
println(n)

# Remove zero rows and columns
solver.A, solver.b, rows_idx, cols_idx = remove_zero_rows_cols!(solver.A, solver.b)

# Solve the reduced system
@time solver.x = solver.A \ solver.b

# Reconstruct the full solution vector
x_full = zeros(2n)
x_full[cols_idx] = solver.x

# Extract Tw (bulk values)
Tw = x_full[1:n]

# Reshape and plot the solution
x_matrix = reshape(Tw, nx+1, ny+1)

fig = Figure(resolution = (800, 800))
ax = Axis(fig[1, 1]; xlabel = "x", ylabel = "y")
heatmap!(ax, x_matrix, colormap = :viridis)
display(fig)

# Solve the problem
#Fliwer.solve!(solver, Fluide; method=IterativeSolvers.bicgstabl, verbose=false, reltol=1e-20)

# Plot the solution
#plot_solution(solver, mesh, circle, capacity)

# Write the solution to a VTK file
#write_vtk("poisson_2d", mesh, solver)

# Fonction pour trouver la valeur maximale d'un tableau
println(maximum(solver.x))