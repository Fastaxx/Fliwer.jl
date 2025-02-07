using Fliwer
using IterativeSolvers
using LinearAlgebra
using SparseArrays
using SpecialFunctions, LsqFit

### 1D Test Case : Monophasic Unsteady Diffusion Equation inside a moving body
# Define the spatial mesh
nx = 320
lx = 10.
x0 = 0.
domain = ((x0, lx),)
mesh = CartesianMesh((nx,), (lx,), (x0,))

# Define the time mesh
# Example of adaptive time stepping
function create_adaptive_times(Tend::Float64, nsteps::Int, ratio::Float64=1.1)
    times = zeros(Float64, nsteps+1)
    times[1] = 0.0
    dt = Tend * (1 - ratio) / (1 - ratio^nsteps)
    for i in 2:nsteps+1
        times[i] = times[i-1] + dt * ratio^(i-2)
    end
    return times
end

Δt = 0.5 * (lx/nx)^2
Tend = 1.0
nt = Int(Tend/Δt)
t = [i*Δt for i in 0:nt]

# Define the body
xf = 0.11*lx   # Interface position
c = 3.0     # Interface velocity
initial_body = Body((x,_=0)->(x - xf), (x,_)->(x), domain, false)  # Initial body
body = Body((x,t, _=0)->(x - xf - c*sqrt(t)), (x,)->(x,), domain, false)  # Body moving to the right
final_body = Body((x,_=0)->(x - xf - c*sqrt(Tend)), (x,_)->(x), domain, false)  # Final body

# Define the space-time mesh
spaceTimeMesh = CartesianSpaceTimeMesh(mesh, t[1:2])

# Identify cells
identify!(mesh, initial_body)
spaceTimeMesh.tag = mesh.tag

# Define the capacity
capacity = Capacity(body, spaceTimeMesh)
capacity_init = Capacity(initial_body, mesh)
capacity_final = Capacity(final_body, mesh)

# Define the operators
operator = SpaceTimeOps(capacity.A, capacity.B, capacity.V, capacity.W, (nx+1, 2))
operator_init = DiffusionOps(capacity_init.A, capacity_init.B, capacity_init.V, capacity_init.W, (nx+1,))
operator_final = DiffusionOps(capacity_final.A, capacity_final.B, capacity_final.V, capacity_final.W, (nx+1,))

# Define the boundary conditions
bc = Dirichlet(0.0)
bc1 = Dirichlet(0.0)

bc_b = BorderConditions(Dict{Symbol, AbstractBoundary}(:top => Dirichlet(0.0), :bottom => Dirichlet(1.0)))

# Define the source term
f = (x,y,z,t)-> 0.0 #sin(x)*cos(10*y)

Fluide = Phase(capacity, operator, f, 1.0)

# Initial condition
u0ₒ = zeros((nx+1))
u0ᵧ = zeros((nx+1))
u0 = vcat(u0ₒ, u0ᵧ)

# Define the solver
solver = MovingDiffusionUnsteadyMono(Fluide, bc_b, bc, Δt, Tend, u0, "CN")

# Solve the problem
solve_MovingDiffusionUnsteadyMono!(solver, Fluide, u0, Δt, Tend, nt, bc_b, bc, body, mesh, t, "CN"; method=Base.:\)

# Write the solution to a VTK file
#write_vtk("moving_heat_1d", mesh, solver)

# Plot the solution
plot_solution(solver, mesh, body, capacity; state_i=1)

# Animation
#animate_solution(solver, mesh, body)

# Analytical solution
function stefan_1d_1ph_analytical(x::Float64)
    t = Tend
    λ = c/2
    return 1.0 - 1.0/erf(2λ) * erf(x/(2*sqrt(t)))
end

function grad_stefan_1d_1ph_analytical(x::Float64)
    t = Tend
    λ = c/2
    return -1.0/(sqrt(t)) * exp(-x^2/(4*t)) / (sqrt(pi) * erf(2λ))
end

using CairoMakie

x=range(x0, stop=lx, length=nx+1)
ls = [final_body.sdf(x[i]) for i in 1:nx+1]
xfaces = x[1:end-1] .+ 0.5*diff(x)
y=[stefan_1d_1ph_analytical(x[i]) for i in 1:nx+1]
y_p=[grad_stefan_1d_1ph_analytical(xfaces[i]) for i in 1:nx]
y[x .> xf + c*Tend] .= 0.0

∇_num = ∇(operator_final, solver.x)

fig = Figure()
ax = Axis(fig[1, 1], xlabel = "x", ylabel = "u", title = "1D 1 phase Stefan problem")
lines!(ax, x, y, color = :blue, linewidth = 2, label = "Analytical solution")
scatter!(ax, x, solver.states[end][1:nx+1], color = :red, label = "Numerical solution")
# add a vertical line to show the interface position
vlines!(ax, xf + c*Tend, color = :black, linestyle = :dash, label = "Interface position")
axislegend(ax)
display(fig)

#u_ana, u_num, global_err, full_err, cut_err, empty_err = check_convergence(stefan_1d_1ph_analytical, solver, capacity_init, mesh, 2, false)

# Plot gradient
fig = Figure()
ax = Axis(fig[1, 1], xlabel = "x", ylabel = "∇u", title = "1D 1 phase Stefan problem - Gradient")
lines!(ax, xfaces, y_p, color = :blue, linewidth = 2, label = "Analytical gradient")
scatter!(ax, x, ∇_num, color = :red, label = "Numerical gradient")
# add a vertical line to show the interface position
vlines!(ax, xf + c*Tend, color = :black, linestyle = :dash, label = "Interface position")
axislegend(ax)
display(fig)

readline()
# Get the last-non zero value of ∇_num
last_non_zero = findlast(x->x!=0.0, ∇_num)
println("Last non zero value of ∇_num: ", ∇_num[last_non_zero])

# Plot Convergence of the L2 error Nusselt number
nx = [20, 40, 80, 160, 320]
dx = 1.0 ./ nx
pseudo_nusselt = [-0.03769354742054475, -0.022132766314183348, -0.017994590389922946, -0.015462342312831003, -0.014356966070581472]
L2_err = [abs(pseudo_nusselt[i] - pseudo_nusselt[end])/abs(pseudo_nusselt[end]) for i in 1:length(pseudo_nusselt)]

fig = Figure()
ax = Axis(fig[1, 1], xlabel = "h", ylabel = "L2 error", title = "Convergence of the L2 error Pseudo Nusselt number - Log scale")
scatter!(ax, log10.(dx), log10.(L2_err), color = :blue, label = "Pseudo Nusselt number (1.43)")
lines!(ax, log10.(dx), log10.(L2_err), color = :black)
lines!(ax, log10.(dx), log10.(dx.^2), color = :black, linestyle = :dash, label = "O(h²)")
lines!(ax, log10.(dx), log10.(dx.^1), color = :black, linestyle = :dashdot, label = "O(h)")
axislegend(ax, position =:rb)
display(fig)

# Model for curve_fit
function fit_model(x, p)
    p[1] .* x .+ p[2]
end

# Fit each on log scale: log(err) = p*log(h) + c
log_h = log.(dx[1:end-1])

function do_fit(log_err)
    fit_result = curve_fit(fit_model, log_h, log_err, [-1.0, 0.0])
    return fit_result.param[1], fit_result.param[2]  # (p_est, c_est)
end

p_global, _ = do_fit(log.(L2_err[1:end-1]))

# Round
p_global = round(p_global, digits=2)

println("Estimated order of convergence (global) = ", p_global)


"""
nx = [20, 40, 80, 160, 320]
dx = 1.0 ./ nx
L2_err_all = [0.2861309825029104, 0.1517310842473155, 0.08706086415610075, 0.055717806923711165, 0.0405692827948743]
L2_err_cut = [0.10094309666605161, 0.036899713160507636, 0.017285819564291833, 0.00948225747532321, 0.005782313993268112]
L2_err_full = [0.2677338797827966, 0.14717585771973002, 0.08532757180183773, 0.054905016178484395, 0.04015509371641087]

fig = Figure()
ax = Axis(fig[1, 1], xlabel = "dx", ylabel = "L2 error", title = "Convergence of the L2 error - Log scale")
scatter!(ax, log10.(dx), log10.(L2_err_all), color = :blue, label = "All cells")
scatter!(ax, log10.(dx), log10.(L2_err_cut), color = :red, label = "Cut cells")
scatter!(ax, log10.(dx), log10.(L2_err_full), color = :green, label = "Full cells")

lines!(ax, log10.(dx), log10.(dx.^2), color = :black, linestyle = :dash, label = "Order 2")
lines!(ax, log10.(dx), log10.(dx.^1), color = :black, linestyle = :dashdot, label = "Order 1")
axislegend(ax, position =:rb)
display(fig)
"""