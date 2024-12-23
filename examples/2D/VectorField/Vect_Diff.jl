using Fliwer
using IterativeSolvers

### 2D Test Case : Monophasic Unsteady Diffusion of a Vector Field inside a Disk
# Define the mesh
nx, ny = 70, 60
lx, ly = 1., 1.
x0, y0 = 0., 0.
domain = ((x0, lx), (y0, ly))
mesh = CartesianMesh((nx, ny), (lx, ly), (x0, y0))

# two new staggered grids are defined, each corresponding to a velocity component, one displaced half a cell width and the other one displaced half a cell height.
mesh_u = CartesianMesh((nx-1, ny), (lx - lx/nx, ly), (x0 + lx/(2*nx), y0))
mesh_v = CartesianMesh((nx, ny-1), (lx, ly - ly/ny), (x0, y0 + ly/(2*ny)))

@show size(mesh_u.centers[1]), size(mesh_u.centers[2])
@show size(mesh_v.centers[1]), size(mesh_v.centers[2])
# Define the body
radius, center = ly/4, (lx/2, ly/2) #.+ (0.01, 0.01)
circle = Body((x,y,_=0)->(sqrt((x-center[1])^2 + (y-center[2])^2) - radius), (x,y,_)->(x,y), domain, false)

# Identify cells
identify!(mesh, circle)
identify!(mesh_u, circle)
identify!(mesh_v, circle)

# Define the capacity
capacity = Capacity(circle, mesh)
capacity_u = Capacity(circle, mesh_u)
capacity_v = Capacity(circle, mesh_v)

# Define the operators
operator_u = DiffusionOps(capacity_u.A, capacity_u.B, capacity_u.V, capacity_u.W, (nx, ny+1))
operator_v = DiffusionOps(capacity_v.A, capacity_v.B, capacity_v.V, capacity_v.W, (nx+1, ny))

# Define the boundary conditions for each velocity component
# For u-component
bc_u = BorderConditions(Dict{Symbol, AbstractBoundary}(:left => Dirichlet(0.0), :right => Dirichlet(0.0), :top => Dirichlet(0.0), :bottom => Dirichlet(0.0)))
ic_u = Dirichlet(1.0)
# For v-component
bc_v = BorderConditions(Dict{Symbol, AbstractBoundary}(:left => Dirichlet(0.0), :right => Dirichlet(0.0), :top => Dirichlet(0.0), :bottom => Dirichlet(0.0)))
ic_v = Dirichlet(0.0)

# Define the source term
fu = (x,y,z,t)-> 2.0
fv = (x,y,z,t)-> 1.0

# Define the phase
Fluide = VectorPhase((capacity_u, capacity_v), (operator_u, operator_v), (fu,fv), 1.0)

# Initial condition
uxₒ0 = zeros(nx*(ny+1))
uyₒ0 = zeros((nx+1)*ny)

uxᵧ0 = zeros(nx*(ny+1))
uyᵧ0 = zeros((nx+1)*ny)

u0 = vcat(uxₒ0, uxᵧ0, uyₒ0, uyᵧ0)

u0x = (uxₒ0, uxᵧ0)
u0y = (uyₒ0, uyᵧ0)

# Define the solver
Δt = 0.01
Tend = 1.0
solver = DiffusionVecUnsteadyMono(Fluide, (bc_u, bc_v), (ic_u, ic_v), Δt, Tend, u0x, u0y)

# Solve the problem
solve_DiffusionVecUnsteadyMono!(solver, Fluide, u0x, u0y, Δt, Tend, (bc_u, bc_v), (ic_u, ic_v))


using CairoMakie

function plot_solution_vector(solver, nx, ny; iter=length(solver.states), title_prefix="")
    # 1) Extract final solution (last iteration)
    sol = solver.states[iter]

    # 2) Decompose solver.x into the four sub-vectors
    #    consistent with the initial condition layout: [uxₒ0; uxᵧ0; uyₒ0; uyᵧ0]
    #    Sizes:
    #      - ux sub-vectors have length nx*(ny+1)
    #      - uy sub-vectors have length (nx+1)*ny
    len_ux = nx*(ny+1)
    len_uy = (nx+1)*ny

    uxₒ  = sol[1 : len_ux]
    uxᵧ  = sol[len_ux + 1 : 2*len_ux]
    uyₒ  = sol[2*len_ux + 1 : 2*len_ux + len_uy]
    uyᵧ  = sol[2*len_ux + len_uy + 1 : end]

    # 3) Reshape them as 2D arrays for plotting
    #    - For ux: shape is (nx, ny+1)
    #    - For uy: shape is (nx+1, ny)
    Uxₒ  = reshape(uxₒ, nx,  ny+1)
    Uxᵧ  = reshape(uxᵧ, nx,  ny+1)
    Uyₒ  = reshape(uyₒ, nx+1, ny)
    Uyᵧ  = reshape(uyᵧ, nx+1, ny)

    # 4) Create a figure with four subplots
    fig = Figure(resolution = (1200, 800))

    # Plot Uxₒ
    ax1 = Axis(fig[1, 1], title = "$title_prefix : Uxₒ")
    hm1 = heatmap!(ax1, Uxₒ, colormap=:viridis)
    Colorbar(fig[1, 2], hm1, label="Uxₒ")

    # Plot Uxᵧ
    ax2 = Axis(fig[1, 3], title = "$title_prefix : Uxᵧ")
    hm2 = heatmap!(ax2, Uxᵧ, colormap=:viridis)
    Colorbar(fig[1, 4], hm2, label="Uxᵧ")

    # Plot Uyₒ
    ax3 = Axis(fig[2, 1], title = "$title_prefix : Uyₒ")
    hm3 = heatmap!(ax3, Uyₒ, colormap=:viridis)
    Colorbar(fig[2, 2], hm3, label="Uyₒ")

    # Plot Uyᵧ
    ax4 = Axis(fig[2, 3], title = "$title_prefix : Uyᵧ")
    hm4 = heatmap!(ax4, Uyᵧ, colormap=:viridis)
    Colorbar(fig[2, 4], hm4, label="Uyᵧ")

    display(fig)
end

# Example usage after solving:
plot_solution_vector(solver, nx, ny, iter=19,title_prefix="DiffusionVecUnsteadyMono")