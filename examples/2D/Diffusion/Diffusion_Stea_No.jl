using Fliwer
using IterativeSolvers

### 2D Test Case : Monophasic Unsteady Diffusion Equation with a NoBody struct
# Define the mesh
nx, ny = 80, 80
lx, ly = 4., 4.
x0, y0 = 0., 0.
domain = ((x0, lx), (y0, ly))
mesh = CartesianMesh((nx, ny), (lx, ly), (x0, y0))

# Define the body
circle = NoBody2D(domain)

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
f = (x,y,z)->1.0

# Define the phase
Fluide = Phase(capacity, operator, f, 1.0)

# Define the solver
solver = DiffusionSteadyMono(Fluide, bc_b, bc)

# Solve the problem
solve_DiffusionSteadyMono!(solver, Fluide; method=Base.:\)

# Plot the solution
#plot_solution(solver, mesh, circle, capacity; state_i=10)

# Animation
#animate_solution(solver, mesh, circle)

# Analytical solution
function u_analytical(x,y)
    # Sum from m,n=1 to 100 m,n odd of 16 L²/((π^4 mn) * (m² + n²)) * sin(mπx/L) * sin(nπy/L)
    sum = 0.0
    for m in 1:2:500
        for n in 1:2:500
            sum += 16 * lx^2 / (π^4 * m * n * (m^2 + n^2)) * sin(m*π*x/lx) * sin(n*π*y/ly)
        end
    end
    return sum
end

u_ana, u_num, global_err, full_err, cut_err, empty_err = check_convergence(u_analytical, solver, capacity, 2)

u_ana[capacity.cell_types .== 0] .= NaN
u_ana = reshape(u_ana, (nx+1, ny+1))

u_num[capacity.cell_types .== 0] .= NaN
u_num = reshape(u_num, (nx+1, ny+1))

err = u_ana - u_num

using CairoMakie
fig = Figure()
ax1 = Axis(fig[1, 1], xlabel = "x", ylabel="y", title="Analytical solution")
ax2 = Axis(fig[1, 2], xlabel = "x", ylabel="y", title="Numerical solution")
hm1 = heatmap!(ax1, u_ana, colormap=:viridis)
hm2 = heatmap!(ax2, u_num, colormap=:viridis)
Colorbar(fig[1, 3], hm1, label="u(x)")
Colorbar(fig[1, 4], hm2, label="u_num(x)")
display(fig)
readline()

# Plot error heatmap
err = reshape(err, (nx+1, ny+1))
fig = Figure()
ax = Axis(fig[1, 1], xlabel = "x", ylabel="y", title="Log error")
hm = heatmap!(ax, log10.(abs.(err)), colormap=:viridis)
Colorbar(fig[1, 2], hm, label="log10(|u(x) - u_num(x)|)")
display(fig)
