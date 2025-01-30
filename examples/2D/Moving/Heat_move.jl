using Fliwer
using IterativeSolvers, SpecialFunctions
using Roots


### 2D Test Case : Monophasic Unsteady Diffusion Equation inside a moving Disk
# Define the mesh
nx, ny = 80, 80
lx, ly = 4., 4.
x0, y0 = 0., 0.
domain = ((x0, lx), (y0, ly))
mesh = CartesianMesh((nx, ny), (lx, ly), (x0, y0))

# Define the time mesh
Δt = 0.25 * (lx/nx)^2
Tend = 0.01
nt = Int(round(Tend/Δt))
t = [i*Δt for i in 0:nt]

# Define the body : Translating disk : radius = 0.5, center = (lx/2 + 0.1*t, ly/2)
radius, center = ly/4, (lx/2, ly/2) .+ (0.01, 0.01)
c = 0.0
initial_body = Body((x, y,_=0) -> sqrt((x - center[1])^2 + (y - center[2])^2) - radius, (x,y,_)->(x,y), domain, false)
body = Body((x, y, t) -> sqrt((x - center[1] - c*t)^2 + (y - center[2])^2) - radius, (x,y,_)->(x,y), domain, false)
final_body = Body((x, y,_=0) -> sqrt((x - center[1] - c*Tend)^2 + (y - center[2])^2) - radius, (x,y,_)->(x,y), domain, false)

# Define the space-time mesh
spaceTimeMesh = CartesianSpaceTimeMesh(mesh, t[1:2])

# Identify cells
identify!(mesh, initial_body)
spaceTimeMesh.tag = mesh.tag

# Define the capacity
capacity = Capacity(body, spaceTimeMesh)
capacity_init = Capacity(initial_body, mesh)
capacity.cell_types = capacity_init.cell_types

capacity_final = Capacity(final_body, mesh)

# CFL log
Vn_1 = capacity.A[2][1:end÷2, 1:end÷2]
Vn = capacity.A[2][end÷2+1:end, end÷2+1:end]
cflmax = 0.2 * Δt./(maximum(Vn_1))
cflmin = 0.2 * Δt./(minimum(Vn_1))

println("CFL number max : ", cflmax)
println("CFL number min : ", cflmin)   

# Define the operators
operator = SpaceTimeOps(capacity.A, capacity.B, capacity.V, capacity.W, (nx+1, ny+1, 2))
operator_init = DiffusionOps(capacity_init.A, capacity_init.B, capacity_init.V, capacity_init.W, (nx+1, ny+1))
operator_final = DiffusionOps(capacity_final.A, capacity_final.B, capacity_final.V, capacity_final.W, (nx+1, ny+1))

# Define the boundary conditions
bc = Dirichlet(0.0)
bc1 = Dirichlet(1.0)

bc_b = BorderConditions(Dict{Symbol, AbstractBoundary}(:left => bc, :right => bc, :top => bc, :bottom => bc))

# Define the source term
f = (x,y,z,t)-> 0.0 #sin(x)*cos(10*y)

Fluide = Phase(capacity, operator, f, 1.0)

# Initial condition
u0ₒ = zeros((nx+1)*(ny+1))
u0ᵧ = ones((nx+1)*(ny+1))
u0 = vcat(u0ₒ, u0ᵧ)

# Define the solver
solver = MovingDiffusionUnsteadyMono2(Fluide, bc_b, bc1, Δt, Tend, u0, "BE")

# Solve the problem
solve_MovingDiffusionUnsteadyMono2!(solver, Fluide, u0, Δt, Tend, nt, bc_b, bc1, body, mesh, t, "BE"; method=Base.:\)

# Plot the solution
plot_solution(solver, mesh, body, capacity; state_i=1)

readline()
# Animation
#animate_solution(solver, mesh, body)

# Analytical solution
function radial_heat_xy(x, y)
    t=nt*Δt
    R=1.0

    function j0_zeros(N; guess_shift=0.25)
        zs = zeros(Float64, N)
        # The m-th zero of J₀ is *roughly* near (m - guess_shift)*π for large m.
        # We'll bracket around that approximate location and refine via find_zero.
        for m in 1:N
            # approximate location
            x_left  = (m - guess_shift - 0.5)*pi
            x_right = (m - guess_shift + 0.5)*pi
            # ensure left>0
            x_left = max(x_left, 1e-6)
            
            # We'll use bisection or Brent's method from Roots.jl
            αm = find_zero(besselj0, (x_left, x_right))
            zs[m] = αm
        end
        return zs
    end

    alphas = j0_zeros(100)
    N=length(alphas)
    r = sqrt((x - center[1] - c*t)^2 + (y - center[2])^2)
    if r >= R
        # Not physically in the domain, so return NaN or handle as you wish.
        return NaN
    end
    
    # If in the disk: sum the series
    s = 0.0
    for m in 1:N
        αm = alphas[m]
        s += exp(-αm^2 * t) * besselj0(αm * (r / R)) / (αm * besselj1(αm))
    end
    return 1.0 - 2.0*s
end

u_ana, u_num, global_err, full_err, cut_err, empty_err = check_convergence(radial_heat_xy, solver, capacity_final, 2)

u_ana[capacity_final.cell_types .== 0] .= NaN
u_ana = reshape(u_ana, (nx+1, ny+1))

u_num[capacity_final.cell_types .== 0] .= NaN
u_num = reshape(u_num, (nx+1, ny+1))

err = u_ana - u_num

using CairoMakie
fig = Figure()
ax1 = Axis(fig[1, 1], xlabel = "x", ylabel="y", title="Analytical solution")
ax2 = Axis(fig[1, 2], xlabel = "x", ylabel="y", title="Numerical solution")
hm1 = heatmap!(ax1, u_ana, colormap=:viridis)
hm2 = heatmap!(ax2, u_num, colormap=:viridis)
Colorbar(fig[1, 3], hm1, label="u(x)")
Colorbar(fig[1, 4], hm2, label="u(x)")
display(fig)
readline()

# Plot error heatmap
err = reshape(err, (nx+1, ny+1))
fig = Figure()
ax = Axis(fig[1, 1], xlabel = "x", ylabel="y", title="Log error")
hm = heatmap!(ax, log10.(abs.(err)), colormap=:viridis)
Colorbar(fig[1, 2], hm, label="log10(|u(x) - u_num(x)|)")
display(fig)
