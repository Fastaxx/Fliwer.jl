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
xf = 0.2*lx   # Interface position
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
u0ₒ = ones((nx+1))
u0ᵧ = zeros((nx+1))
u0 = vcat(u0ₒ, u0ᵧ)

@show u0
# Define the solver
s = MotionDiffusionUnsteadyMono(Fluide, bc_b, bc, Δt, Tend, u0, "BE")
#solve_MotionDiffusionUnsteadyMono!(solver, Fluide, u0, Δt, Tend, nt, bc_b, bc, body, mesh, t, "CN")

# Set iteration parameters
max_iter = 200
tol      = 1e-6

old_xf = xf
iter   = 0
err    = Inf

function run_stefan_iteration!(
    s,
    body::Body,
    mesh::CartesianMesh,
    spaceTimeMesh::CartesianSpaceTimeMesh,
    capacity::Capacity,
    operator::SpaceTimeOps,
    Fluide::Phase,
    bc_b::BorderConditions,
    bc::AbstractBoundary,
    u0::Vector{Float64},
    domain,
    ρ::Float64,
    L::Float64,
    max_iter::Int,
    tol::Float64,
    Δt::Float64,
    Tend::Float64,
    xf::Float64
)
    err::Float64 = Inf
    old_xf::Float64 = xf
    iter::Int = 0

    while (iter < max_iter) && (err > tol)
        iter += 1

        # 1) Solve the linear system
        A_reduced, b_reduced, rows_idx, cols_idx = remove_zero_rows_cols!(s.A, s.b)
        x_reduced = A_reduced \ b_reduced
        s.x = zeros(size(s.A, 1))
        s.x[cols_idx] .= x_reduced
        Tᵢ = s.x

        @show Tᵢ
        # 2) Update volumes / compute new interface
        Vn_1 = capacity.A[2][1:end÷2, 1:end÷2]
        Vn   = capacity.A[2][end÷2+1:end, end÷2+1:end]
        Hₙ   = sum(diag(Vn))
        Hₙ₊₁ = sum(diag(Vn_1))

        # Compute flux
        nx = Fluide.operator.size[1]
        W! = operator.Wꜝ[1:nx, 1:nx]
        G = operator.G[1:nx, 1:nx]
        H = operator.H[1:nx, 1:nx]
        V = operator.V[1:nx, 1:nx]
        Tₒ, Tᵧ = Tᵢ[1:nx], Tᵢ[nx+1:end]
        Interface_term = H' * W! * G * Tₒ + H' * W! * H * Tᵧ
        Interface_term = 1/(ρ*L) * sum(Interface_term)

        # New interface position
        xf_new = Hₙ + Interface_term
        err = abs(xf_new - old_xf)
        println("Iteration $iter | xf = $xf_new | error = $err")

        # 3) Update geometry if not converged
        if err <= tol
            break
        end
        old_xf = xf_new
        xf = xf_new
        ΔH = xf - Hₙ

        # Rebuild domain
        body = Body((xx,t, _=0)->(xx - xf), (xx,)->(xx,), domain, false)
        #identify!(mesh, body)
        spaceTimeMesh.tag = mesh.tag

        capacity = Capacity(body, spaceTimeMesh)
        operator = SpaceTimeOps(capacity.A, capacity.B, capacity.V, capacity.W, (nx, 2))
        Fluide   = Phase(capacity, operator, (x,y,z,t)->0.0, 1.0)

        s = MotionDiffusionUnsteadyMono(Fluide, bc_b, bc, Δt, Tend, u0, "BE")
    end

    if err <= tol
        println("Converged after $iter iterations with xf = $xf, error = $err")
    else
        println("Reached max_iter = $max_iter with xf = $xf, error = $err")
    end

    return xf
end

xf_final = run_stefan_iteration!(
     s, body, mesh, spaceTimeMesh, capacity,
     operator, Fluide, bc_b, bc, u0, domain,
     ρ, L, max_iter, tol, Δt, Tend, xf )

"""
Vin = [Vn[i, i] for i in 1:size(Vn, 1)]
Vin1 = [Vn_1[i, i] for i in 1:size(Vn_1, 1)]

fig = Figure(size = (800, 400))
ax = Axis(fig[1, 1], xlabel = "x", ylabel = "V", title = "Vn and Vn+1")
scatter!(ax, x, Vin, color = :blue, label = "Vn")
scatter!(ax, x, Vin1, color = :red, label = "Vn+1")
axislegend(ax, position = :rt)
display(fig)
"""