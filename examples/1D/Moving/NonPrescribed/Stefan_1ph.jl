using Fliwer
using IterativeSolvers
using LinearAlgebra
using SparseArrays
using SpecialFunctions, LsqFit
using CairoMakie

### 1D Test Case : Monophasic Unsteady Diffusion Equation inside a moving body
# Define the spatial mesh
nx = 160
lx = 10.
x0 = 0.
domain = ((x0, lx),)
mesh = CartesianMesh((nx,), (lx,), (x0,))
x = range(x0, stop = lx, length = nx+1)

Δt = 0.5 * (lx/nx)^2
Tend = 1.0
nt = 3
t = [i*Δt for i in 0:nt]

# Define the body
xf = 0.02*lx   # Interface position
@show xf
initial_body = Body((x,_=0)->(x - xf), (x,_)->(x), domain, false)  # Initial body
body = Body((x,t, _=0)->(x - xf), (x,)->(x,), domain, false)

# Define the space-time mesh
spaceTimeMesh = CartesianSpaceTimeMesh(mesh, t[1:2])

# Identify cells
identify!(mesh, initial_body)
spaceTimeMesh.tag = mesh.tag

# Define the capacity
capacity = Capacity(body, spaceTimeMesh)

# Define the operators
operator = SpaceTimeOps(capacity.A, capacity.B, capacity.V, capacity.W, (nx+1, 2))

# Define the boundary conditions
bc = Dirichlet(0.0)

bc_b = BorderConditions(Dict{Symbol, AbstractBoundary}(:top => Dirichlet(0.0), :bottom => Dirichlet(1.0)))

# Define the source term
f = (x,y,z,t)-> 0.0 

ρ = 1.0
L = 1.0
Fluide = Phase(capacity, operator, f, 1.0)

# Initial condition
u0ₒ = zeros((nx+1))
u0ᵧ = zeros((nx+1))
u0 = vcat(u0ₒ, u0ᵧ)

# Define the solver
s = MotionDiffusionUnsteadyMono(Fluide, bc_b, bc, Δt, Tend, u0, "BE")
solve_MotionDiffusionUnsteadyMono!(s, Fluide, u0, Δt, Tend, nt, bc_b, bc, body, mesh, t, "BE", xf)
