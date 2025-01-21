using Fliwer
using IterativeSolvers
using SparseArrays
using LinearAlgebra
using CairoMakie

### 2D Test Case : Monophasic Steady Diffusion Equation inside a Disk
# Define the mesh
nx, ny = 40, 40
lx, ly = 4., 4.
x0, y0 = 0., 0.
domain = ((x0, lx), (y0, ly))
mesh = CartesianMesh((nx, ny), (lx, ly), (x0, y0))

# Define the body
radius, center = ly/4, (lx/2, ly/2) .+ (0.01, 0.01)
circle = Body((x,y,_=0)->(sqrt((x-center[1])^2 + (y-center[2])^2) - radius), (x,y,_)->(x,y), domain, false)

# Identify cells
identify!(mesh, circle)

# Define the capacity
capacity = Capacity(circle, mesh)

# Define the operators
operator = DiffusionOps(capacity.A, capacity.B, capacity.V, capacity.W, (nx+1, ny+1))

# Define the boundary conditions 
bc = Dirichlet(0.0)
bc1 = Dirichlet(0.0)
bc_neumann = Neumann(1.0)
bc_periodic = Periodic()

bc_b = BorderConditions(Dict{Symbol, AbstractBoundary}(:left => bc1, :right => bc1, :top => bc1, :bottom => bc1))

# Define the source term
f = (x,y,_)-> 4.0 #sin(x)*cos(10*y)
K = (x,y,_) -> begin
    r = sqrt((x - lx/2)^2 + (y - lx/3)^2)
    if r <= 1.0
        return 10.0
    else
        return 0.1
    end
end
K=1.0

Fluide = Phase(capacity, operator, f, K)

# Define the solver
solver = DiffusionSteadyMono(Fluide, bc_b, bc)

# Solve the problem
Fliwer.solve_DiffusionSteadyMono!(solver, Fluide; method=Base.:\)
#Fliwer.solve_DiffusionSteadyMono!(solver, Fluide; method=IterativeSolvers.bicgstabl, verbose=false, reltol=1e-20)

# Plot the solution
plot_solution(solver, mesh, circle, capacity)

# Write the solution to a VTK file
#write_vtk("poisson_2d", mesh, solver)

# Fonction pour trouver la valeur maximale d'un tableau
println(maximum(solver.x))

# Analytical solution
u_analytical(x,y) = 1.0 - (x-center[1])^2 - (y-center[2])^2

u_ana, u_num, global_err, full_err, cut_err, empty_err = check_convergence(u_analytical, solver, capacity, 2)

u_ana[capacity.cell_types .== 0] .= NaN
u_ana = reshape(u_ana, (nx+1, ny+1))

u_num[capacity.cell_types .== 0] .= NaN
u_num = reshape(u_num, (nx+1, ny+1))

err = u_ana .- u_num

using CairoMakie
fig = Figure()
ax1 = Axis(fig[1, 1], xlabel = "x", ylabel="y", title="Analytical solution")
ax2 = Axis(fig[1, 2], xlabel = "x", ylabel="y", title="Numerical solution")
heatmap!(ax1, u_ana, colormap=:viridis)
heatmap!(ax2, u_num, colormap=:viridis)
Colorbar(fig[1, 3], label="u(x)")
display(fig)

# Plot error heatmap
err = reshape(err, (nx+1, ny+1))
fig = Figure()
ax = Axis(fig[1, 1], xlabel = "x", ylabel="y", title="Log error")
hm = heatmap!(ax, log10.(abs.(err)), colormap=:viridis)
Colorbar(fig[1, 2], hm, label="log10(|u(x) - u_num(x)|)")
display(fig)
