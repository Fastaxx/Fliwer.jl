using Fliwer
using IterativeSolvers

### 2D Test Case : Monophasic Unsteady Convection of a Vector Field inside a Disk
# Define the mesh
nx, ny = 80, 80
lx, ly = 4., 4.
x0, y0 = 0., 0.
domain = ((x0, lx), (y0, ly))
mesh = CartesianMesh((nx, ny), (lx, ly), (x0, y0))

# two new staggered grids are defined, each corresponding to a velocity component, one displaced half a cell width and the other one displaced half a cell height.
mesh_u = CartesianMesh((nx-1, ny), (lx - lx/nx, ly), (x0 + lx/(2*nx), y0))
mesh_v = CartesianMesh((nx, ny-1), (lx, ly - ly/ny), (x0, y0 + ly/(2*ny)))

@show size(mesh_u.centers[1]), size(mesh_u.centers[2])
@show mesh_u.centers
@show mesh_v.centers
# Define the body
radius, center = ly/4, (lx/2, ly/2) #.+ (0.01, 0.01)
circle = Body((x,y,_=0)->(sqrt((x-center[1])^2 + (y-center[2])^2) - radius), (x,y,_)->(x,y), domain, false)

# Initial condition
uxₒ0 = zeros(nx*(ny+1))
uyₒ0 = zeros((nx+1)*ny)

uxᵧ0 = zeros(nx*(ny+1))
uyᵧ0 = zeros((nx+1)*ny)

# Fill uxₒ0 using mesh_u and nan values for outside the circle
for j in 1:length(mesh_u.centers[2])
    for i in 1:length(mesh_u.centers[1])
        x = mesh_u.centers[1][i]
        y = mesh_u.centers[2][j]
        idx = i + (j - 1) * (length(mesh_u.centers[1])+1)
        uxₒ0[idx] = -0.5 * (y - center[2])   # rigid rotation around the origin
        if sqrt((mesh_u.centers[1][i] - center[1])^2 + (mesh_u.centers[2][j] - center[2])^2) > radius
            idx = i + (j - 1) * (length(mesh_u.centers[1])+1)
            #uxₒ0[idx] = NaN
        end
    end
end

# Fill uyₒ0 using mesh_v
for j in 1:length(mesh_v.centers[2])
    for i in 1:length(mesh_v.centers[1])
        x = mesh_v.centers[1][i]
        y = mesh_v.centers[2][j]
        idx = i + (j - 1) * (length(mesh_v.centers[1])+1)
        uyₒ0[idx] = 0.5 * (x - center[1])   # rigid rotation around the origin
        if sqrt((mesh_v.centers[1][i] - center[1])^2 + (mesh_v.centers[2][j] - center[2])^2) > radius
            idx = i + (j - 1) * (length(mesh_v.centers[1])+1)
            #uyₒ0[idx] = NaN
        end
    end
end

"""
# Plot the initial condition
using CairoMakie

fig = Figure(size = (800, 800))
ax1 = Axis(fig[1, 1], title = "uxₒ0", aspect=DataAspect())
hm1 = heatmap!(ax1, reshape(uxₒ0, nx, ny+1), colormap=:viridis)
Colorbar(fig[1, 2], hm1, label="uxₒ0")

ax2 = Axis(fig[2, 1], title = "uyₒ0", aspect=DataAspect())
hm2 = heatmap!(ax2, reshape(uyₒ0, nx+1, ny), colormap=:viridis)
Colorbar(fig[2, 2], hm2, label="uyₒ0")

display(fig)
readline()
"""

u0 = vcat(uxₒ0, uxᵧ0, uyₒ0, uyᵧ0)

u0x = (uxₒ0, uxᵧ0)
u0y = (uyₒ0, uyᵧ0)

# Identify cells
identify!(mesh, circle)
identify!(mesh_u, circle)
identify!(mesh_v, circle)

# Define the capacity
capacity = Capacity(circle, mesh)
capacity_u = Capacity(circle, mesh_u)
capacity_v = Capacity(circle, mesh_v)

# Define the operators
#operator = AdvectionVecOps((capacity_u, capacity_v), (nx-1,ny-1), (uxₒ0, uyₒ0), (uxᵧ0, uyᵧ0))
operator_u = ConvectionOps(capacity_u.A, capacity_u.B, capacity_u.V, capacity_u.W, (nx, ny+1), (uxₒ0, uyₒ0), vcat(uxᵧ0, uxᵧ0))
operator_v = ConvectionOps(capacity_v.A, capacity_v.B, capacity_v.V, capacity_v.W, (nx+1, ny), (uyₒ0, uyₒ0), vcat(uyᵧ0, uyᵧ0))

# Define the boundary conditions for each velocity component
# For u-component
bc_u = BorderConditions(Dict{Symbol, AbstractBoundary}(:left => Dirichlet(0.0), :right => Dirichlet(0.0), :top => Dirichlet(0.0), :bottom => Dirichlet(0.0)))
ic_u = Dirichlet(1.0)

# For v-component
bc_v = BorderConditions(Dict{Symbol, AbstractBoundary}(:left => Dirichlet(0.0), :right => Dirichlet(0.0), :top => Dirichlet(0.0), :bottom => Dirichlet(0.0)))
ic_v = Dirichlet(0.0)

# Define the source term
fu = (x,y,z,t)-> begin
    r = sqrt((x - lx/2)^2 + (y - lx/4)^2)
    if r <= 0.4
        return 1.0
    else
        return 0.0
    end
    end
fv = (x,y,z,t)-> 0.0

# Define the phase
Fluide = VectorPhase((capacity_u, capacity_v), (operator_u, operator_v), (fu,fv), 1.0)

# Define the solver
Δt = 0.01
Tend = 1.0
solver = ConvectionVecUnsteadyMono(Fluide, (bc_u, bc_v), (ic_u, ic_v), Δt, Tend, u0)

# Solve the problem
solve_ConvectionVecUnsteadyMono!(solver, Fluide, u0, Δt, Tend, (bc_u, bc_v), (ic_u, ic_v), method=Base.:\)

write_vtk("ConvectionVecUnsteadyMono", mesh, solver)