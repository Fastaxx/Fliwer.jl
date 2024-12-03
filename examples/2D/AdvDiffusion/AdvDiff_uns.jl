using Fliwer
using IterativeSolvers
using LinearAlgebra

### 2D Test Case : Monophasic Unsteady Advection-Diffusion Equation inside a Disk
# Define the mesh
nx, ny = 160, 160
lx, ly = 16., 16.
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

# Initialize the velocity field by a constant value
# Number of nodes
N = (nx + 1) * (ny + 1)

# Initialize velocity components
uₒx = zeros(N)
uₒy = 10.0 .* ones(N) 

# Interleave the components
uₒ = Vector{Float64}(undef, 2N)
uₒ[1:2:2N] .= uₒx
uₒ[2:2:2N] .= uₒy

# For boundary velocities, if they are zero:
uᵧ = zeros(2N)

# Define the operators
operator = ConvectionOps(capacity.A, capacity.B, capacity.V, capacity.W, (nx+1, ny+1), uₒ, uᵧ)

# Define the boundary conditions
ic = Dirichlet(0.0)
bc = Dirichlet(0.0)

bc_b = BorderConditions(Dict{Symbol, AbstractBoundary}(:left => bc, :right => bc, :top => bc, :bottom => bc))

# Define the source term
f = (x,y,z,t)-> 0.0 #sin(x)*cos(10*y)

Fluide = Phase(capacity, operator, f, 1.0)

# Initial condition
# Analytical solution T (x, t) = 10/(4D (t + 1/2)) exp(− ‖x−xc(t)‖^2/(4D(t+ 1/2)
T0ₒ = zeros((nx+1)*(ny+1))
T0ᵧ = zeros((nx+1)*(ny+1))
x, y = mesh.centers[1], mesh.centers[2]

for j in 1:ny
    for i in 1:nx
        T0ₒ[i + (j-1)*(nx+1)] = 10/(4*1/2) * exp(-((x[i]-center[1])^2 + (y[j]-center[2])^2)/(4*1/2))
        T0ᵧ[i + (j-1)*(nx+1)] = 10/(4*1/2) * exp(-((x[i]-center[1])^2 + (y[j]-center[2])^2)/(4*1/2))
    end
end

# Initialize temperature arrays with zeros
T0ₒ = zeros((nx + 1) * (ny + 1))
T0ᵧ = zeros((nx + 1) * (ny + 1))

# Coordinates of the grid centers
x_coords = mesh.centers[1]
y_coords = mesh.centers[2]

# Define the half-width of the square region (in number of cells)
square_half_width = 5  # Adjust this value to control the size of the square

# Find the indices of the center cell
center_i = findfirst(x -> x >= center[1], x_coords)
center_j = findfirst(y -> y >= center[2], y_coords)

# Determine the range of indices around the center
i_min = max(center_i - square_half_width, 1)
i_max = min(center_i + square_half_width, nx + 1)
j_min = max(center_j - square_half_width, 1)
j_max = min(center_j + square_half_width, ny + 1)

# Set T = 1 in the square region at the center
for j in j_min:j_max
    for i in i_min:i_max
        idx = i + (j - 1) * (nx + 1)
        T0ₒ[idx] = 1.0
        T0ᵧ[idx] = 1.0
    end
end

# Combine the temperature arrays if needed
T0 = vcat(T0ₒ, T0ᵧ)

# Initialize temperature arrays with zeros
T0ₒ = zeros((nx + 1) * (ny + 1))
T0ᵧ = zeros((nx + 1) * (ny + 1))

# Coordinates of the grid centers
x_coords = mesh.centers[1]
y_coords = mesh.centers[2]

# Define the radius of the circle (in physical units)
circle_radius = lx / 60  # Adjust this value to control the size of the circle

# Loop over all grid points
for j in 1:(ny )
    for i in 1:(nx )
        idx = i + (j - 1) * (nx + 1)
        x_i = x_coords[i]
        y_j = y_coords[j]
        # Compute the distance from the center
        distance = sqrt((x_i - center[1])^2 + (y_j - center[2])^2)
        # If the point is inside the circle, set T=1
        if distance <= circle_radius
            T0ₒ[idx] = 1.0
            T0ᵧ[idx] = 1.0
        end
    end
end

# Combine the temperature arrays if needed
T0 = vcat(T0ₒ, T0ᵧ)

# Define the solver
Δt = 0.01
Tend = 1.0
solver = AdvectionDiffusionUnsteadyMono(Fluide, bc_b, ic, Δt, Tend, T0)

# Solve the problem
solve_AdvectionDiffusionUnsteadyMono!(solver, Fluide, T0, Δt, Tend, bc_b, ic; method=IterativeSolvers.bicgstabl, verbose=false, reltol=1e-20)

# Plot the solution
plot_solution(solver, mesh, circle, capacity)

# Write the solution to a VTK file
write_vtk("advdiff", mesh, solver)

plot_profile(solver, mesh; x=lx/2.01)

# Animation
#animate_solution(solver, mesh, circle)