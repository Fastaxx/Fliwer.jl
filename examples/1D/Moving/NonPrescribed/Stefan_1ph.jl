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
nt = 2
t = [i*Δt for i in 0:nt]

# Define the body
xf = 0.02*lx   # Interface position
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
solver = MotionDiffusionUnsteadyMono(Fluide, bc_b, bc, Δt, Tend, u0, "CN")
solve_MotionDiffusionUnsteadyMono!(solver, Fluide, u0, Δt, Tend, nt, bc_b, bc, body, mesh, t, "CN")

"""
# Newton initial
Vn_1 = capacity.A[2][1:end÷2, 1:end÷2]
Vn = capacity.A[2][end÷2+1:end, end÷2+1:end]

# Height function Hn and Hn+1,k : Sum of the Volume of the cells
Hⁿ = sum(diag(Vn))
Hⁿ⁺¹ = sum(diag(Vn_1))
ΔH = Hⁿ⁺¹ - Hⁿ

# Interface term : 
W! = operator.Wꜝ[1:nx+1, 1:nx+1]
G = operator.G[1:nx+1, 1:nx+1]
H = operator.H[1:nx+1, 1:nx+1]
V = operator.V[1:nx+1, 1:nx+1]
Interface_term = H' * W! * G * u0ₒ + H' * W! * H * u0ᵧ
Interface_term = -1/(ρ*L) * sum(Interface_term)
@show Interface_term

# Update the height function
res = Hⁿ⁺¹
Hⁿ⁺¹ = Hⁿ + Interface_term

@show res
@show Hⁿ⁺¹


Vin = [Vn[i, i] for i in 1:size(Vn, 1)]
Vin1 = [Vn_1[i, i] for i in 1:size(Vn_1, 1)]

fig = Figure(size = (800, 400))
ax = Axis(fig[1, 1], xlabel = "x", ylabel = "V", title = "Vn and Vn+1")
scatter!(ax, x, Vin, color = :blue, label = "Vn")
scatter!(ax, x, Vin1, color = :red, label = "Vn+1")
axislegend(ax, position = :rt)
display(fig)
"""