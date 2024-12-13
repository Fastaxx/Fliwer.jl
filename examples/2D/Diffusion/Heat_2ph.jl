using Fliwer
using IterativeSolvers
using WriteVTK

### 2D Test Case : Diphasic Unsteady Diffusion Equation with a Disk
# Define the mesh
nx, ny = 120, 120
lx, ly = 8., 8.
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

ic = InterfaceConditions(ScalarJump(1.0, 2.0, 0.0), FluxJump(1.0, 1.0, 0.0))

# Define the source term
f1 = (x,y,z,t)->0.0
f2 = (x,y,z,t)->0.0

# Define the phases
Fluide_1 = Phase(capacity, operator, f1, 1.0)
Fluide_2 = Phase(capacity_c, operator_c, f2, 1.0)

# Initial condition
u0ₒ1 = ones((nx+1)*(ny+1))
u0ᵧ1 = ones((nx+1)*(ny+1))
u0ₒ2 = zeros((nx+1)*(ny+1))
u0ᵧ2 = zeros((nx+1)*(ny+1))
u0 = vcat(u0ₒ1, u0ᵧ1, u0ₒ2, u0ᵧ2)

# Define the solver
Δt = 0.01
Tend = 1.0
solver = DiffusionUnsteadyDiph(Fluide_1, Fluide_2, bc_b, ic, Δt, Tend, u0)

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

# Store the size of the original system
n = Int(size(solver.A, 1) /4)
println(n)

# Remove zero rows and columns
solver.A, solver.b, rows_idx, cols_idx = remove_zero_rows_cols!(solver.A, solver.b)

# Solve the reduced system
@time solver.x = solver.A \ solver.b

# Reconstruct the full solution vector
x_full = zeros(4n)
x_full[cols_idx] = solver.x

# Extract Tw (bulk values)
Tw1 = x_full[1:n]
Tg1 = x_full[n+1:2n]
Tw2 = x_full[2n+1:3n]
Tg2 = x_full[3n+1:4n]

# Reshape and plot the solution
Tw1_matrix = reshape(Tw1, nx+1, ny+1)
Tg1_matrix = reshape(Tg1, nx+1, ny+1)
Tw2_matrix = reshape(Tw2, nx+1, ny+1)
Tg2_matrix = reshape(Tg2, nx+1, ny+1)

fig = Figure()
ax1 = Axis(fig[1, 1]; xlabel = "x", ylabel = "y")
ax2 = Axis(fig[1, 2]; xlabel = "x", ylabel = "y")
ax3 = Axis(fig[2, 1]; xlabel = "x", ylabel = "y")
ax4 = Axis(fig[2, 2]; xlabel = "x", ylabel = "y")
heatmap!(ax1, Tw1_matrix, colormap = :viridis)
heatmap!(ax2, Tg1_matrix, colormap = :viridis)
heatmap!(ax3, Tw2_matrix, colormap = :viridis)
heatmap!(ax4, Tg2_matrix, colormap = :viridis)
display(fig)


# Solve the problem
#solve_DiffusionUnsteadyDiph!(solver, Fluide_1, Fluide_2, u0, Δt, Tend, bc_b, ic; method=IterativeSolvers.gmres, maxiter=10000, verbose=false)

# Write the solution to a VTK file
#write_vtk("solution", mesh, solver)

# Plot the solution
#plot_solution(solver, mesh, circle, capacity)

# Plot the Profile
#plot_profile(solver, mesh; x=lx/2.01)

# Animation
#animate_solution(solver, mesh, circle)
