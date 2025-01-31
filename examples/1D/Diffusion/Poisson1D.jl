using Fliwer
using IterativeSolvers
using LinearAlgebra, SparseArrays

### 1D Test Case : Monophasic Steady Diffusion Equation
# Define the mesh
nx = 640
lx = 1.
x0 = 0.
domain = ((x0, lx),)
mesh = CartesianMesh((nx,), (lx,), (x0,))

# Define the body
center, radius = 0.5, 0.21
body = Body((x,_=0)->sqrt((x - center)^2) - radius, (x,)->(x,), domain, false)

# Identify cells
identify!(mesh, body)

# Define the capacity
capacity = Capacity(body, mesh)


@show capacity.V
# Redefine W and V
pₒ = [capacity.C_ω[i][1] for i in 1:length(capacity.C_ω)]
pᵧ = [capacity.C_γ[i][1] for i in 1:length(capacity.C_ω)]
p = vcat(pₒ, pᵧ)
grad = ∇(operator, p)
W_new = [grad[i] * capacity.W[1][i,i] for i in 1:length(grad)]
W_new = spdiagm(0 => W_new)

pₒ = [(capacity.C_ω[i][1]^2)/2 for i in 1:length(capacity.C_ω)]
pᵧ = [(capacity.C_γ[i][1]^2)/2 for i in 1:length(capacity.C_ω)]

p = vcat(pₒ, pᵧ)
grad = ∇(operator, p)

qω = vcat(grad)
qγ = vcat(grad)

div = ∇_(operator, qω, qγ)
V_new = spdiagm(0 => div)

capacity.W = (W_new,)
capacity.V = V_new

@show capacity.V
# Define the operators
operator = DiffusionOps(capacity.A, capacity.B, capacity.V, capacity.W, (nx+1,))

# Define the boundary conditions
bc = Dirichlet(0.0)
bc1 = Dirichlet(0.0)

bc_b = BorderConditions(Dict{Symbol, AbstractBoundary}(:top => bc1, :bottom => bc1))

# Define the source term
f = (x, y, _=0) -> 1.0

Fluide = Phase(capacity, operator, f, 1.0)

# Define the solver
solver = DiffusionSteadyMono(Fluide, bc_b, bc)

# Solve the problem
solve_DiffusionSteadyMono!(solver, Fluide; method=Base.:\)

# Plot the solution
plot_solution(solver, mesh, body, capacity)

# Analytical solution
a, b = center - radius, center + radius
u_analytical = (x) -> -1/2 * (x-a)*(x-b)

u_ana, u_num, err, global_err, full_err, cut_err, empty_err = check_convergence(u_analytical, solver, capacity, 2)

# Plot the error
err[capacity.cell_types .== 0] .= NaN
using CairoMakie
fig = Figure()
ax = Axis(fig[1, 1], xlabel="x", ylabel="Log10 of the error", title="Monophasic Steady Diffusion Equation")
scatter!(ax, log10.(abs.(err)), label="Log10 of the error")
display(fig)