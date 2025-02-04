using Fliwer
using IterativeSolvers
using SparseArrays
using LinearAlgebra
using CairoMakie

### 2D Test Case : Monophasic Steady Diffusion Equation inside a Disk
# Define the mesh
nx, ny = 160, 160
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

function volume_redefinitio!(capacity::Capacity{2}, operator::DiffusionOps)
    pₒ = [capacity.C_ω[i][1] for i in 1:length(capacity.C_ω)]
    pᵧ = [capacity.C_γ[i][1] for i in 1:length(capacity.C_ω)]
    p = vcat(pₒ, pᵧ)
    grad = ∇(operator, p)
    grad_x, grad_y = grad[1:end÷2], grad[end÷2+1:end]
    W_x_new = [grad_x[i] * capacity.W[1][i,i] for i in 1:length(grad_x)]
    W_y_new = [grad_y[i] * capacity.W[2][i,i] for i in 1:length(grad_y)]
    W_new = (spdiagm(0 => W_x_new), spdiagm(0 => W_y_new))


    pₒx = [(capacity.C_ω[i][1]^2)/2 for i in 1:length(capacity.C_ω)]
    pᵧx = [(capacity.C_γ[i][1]^2)/2 for i in 1:length(capacity.C_ω)]

    pₒy = [(capacity.C_ω[i][2]^2)/2 for i in 1:length(capacity.C_ω)]
    pᵧy = [(capacity.C_γ[i][2]^2)/2 for i in 1:length(capacity.C_ω)]

    px = vcat(pₒx, pᵧx)
    py = vcat(pₒy, pᵧy)
    gradx = ∇(operator, px)
    grady = ∇(operator, py)
    gradx_x, gradx_y = gradx[1:end÷2], gradx[end÷2+1:end]
    grady_x, grady_y = grady[1:end÷2], grady[end÷2+1:end]

    qωx = vcat(gradx_x, gradx_y)
    qγx = vcat(grady_x, grady_y)

    div_x = ∇_(operator, qωx, qγx)
    div_y = ∇_(operator, qωx, qγx)

    @show minimum(div_x), minimum(div_y)
    @show maximum(div_x), maximum(div_y)

    V_new = max.(div_x, div_y)
    V_new = spdiagm(0 => V_new)

    capacity.W = W_new
    capacity.V = V_new
end

volume_redefinitio!(capacity, operator)
operator=DiffusionOps(capacity.A, capacity.B, capacity.V, capacity.W, (nx+1, ny+1))

# Define the boundary conditions 
bc = Dirichlet((x,y,z) -> 0)
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

readline()
# Write the solution to a VTK file
#write_vtk("poisson_2d", mesh, solver)

# Fonction pour trouver la valeur maximale d'un tableau
println(maximum(solver.x))

# Analytical solution
u_analytical(x,y) = 1.0 - (x-center[1])^2 - (y-center[2])^2
∇x_analytical(x,y) = -2*(x-center[1])
∇y_analytical(x,y) = -2*(y-center[2])

u_ana, u_num, global_err, full_err, cut_err, empty_err = check_convergence(u_analytical, solver, capacity, 2, false)

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

readline()

# Gradient error
∇_num = ∇(operator, solver.x)
∇_numx = ∇_num[1:end÷2]
∇_numy = ∇_num[end÷2+1:end]

∇x_ana = [∇x_analytical(x,y) for x in range(x0, stop=lx, length=nx+1), y in range(y0, stop=ly, length=ny+1)]
∇y_ana = [∇y_analytical(x,y) for x in range(x0, stop=lx, length=nx+1), y in range(y0, stop=ly, length=ny+1)]

∇x_ana[capacity.cell_types .== 0] .= NaN
∇x_ana = reshape(∇x_ana, (nx+1, ny+1))

∇y_ana[capacity.cell_types .== 0] .= NaN
∇y_ana = reshape(∇y_ana, (nx+1, ny+1))

∇_numx[capacity.cell_types .== 0] .= NaN
∇_numx = reshape(∇_numx, (nx+1, ny+1))

∇_numy[capacity.cell_types .== 0] .= NaN
∇_numy = reshape(∇_numy, (nx+1, ny+1))

∇x_num = reshape(∇_numx, (nx+1, ny+1))
∇y_num = reshape(∇_numy, (nx+1, ny+1))

∇x_err = ∇x_ana .- ∇x_num
∇y_err = ∇y_ana .- ∇y_num

using CairoMakie
fig = Figure()
ax1 = Axis(fig[1, 1], xlabel = "x", ylabel="y", title="Analytical gradient in x")
ax2 = Axis(fig[1, 2], xlabel = "x", ylabel="y", title="Numerical gradient in x")
ax3 = Axis(fig[1, 3], xlabel = "x", ylabel="y", title="Log10 Error in gradient in x")
hm1 = heatmap!(ax1, ∇x_ana, colormap=:viridis)
hm2 = heatmap!(ax2, ∇x_num, colormap=:viridis)
hm3 = heatmap!(ax3, log10.(abs.(∇x_err)), colormap=:viridis)
Colorbar(fig[1, 4], hm1, label="∇x(x)")
Colorbar(fig[1, 5], hm2, label="∇x(x)")
Colorbar(fig[1, 6], hm3, label="∇x(x)")
display(fig)
