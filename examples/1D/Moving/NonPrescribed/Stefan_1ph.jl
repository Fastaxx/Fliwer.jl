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

function newton_1step!(
    s::Solver,
    phase::Phase,
    Tᵢ::Vector{Float64},
    xf::Float64,
    bc_b::BorderConditions,
    bc::AbstractBoundary,
    mesh::CartesianMesh,
    tspan::NTuple{2,Float64},
    Δt::Float64,
    scheme::String;
    method = IterativeSolvers.gmres,
    max_iter = 1000,
    tol = 1e-15,
    kwargs...
)
    println("Newton routine for one time step...")
    nx, _ = phase.operator.size
    ρ, L = 1.0, 1.0

    residuals = Float64[]
    err = Inf
    old_xf = xf
    iter = 0

    while (iter < max_iter) && (err > tol)
        iter += 1

        # 1) Solve the linear system
        if method == \
            A_reduced, b_reduced, rows_idx, cols_idx = remove_zero_rows_cols!(s.A, s.b)
            x_reduced = A_reduced \ b_reduced
            s.x = zeros(size(s.A, 1))
            s.x[cols_idx] .= x_reduced
        else
            log = get(kwargs, :log, false)
            if log
                s.x, s.ch = method(s.A, s.b; kwargs...)
            else
                s.x = method(s.A, s.b; kwargs...)
            end
        end
        Tᵢ .= s.x

        # 2) Compute new interface position using current phase data
        #    Vn_1: the cell volumes before interface and Vn: after interface
        Vn_1 = phase.capacity.A[2][1:end÷2, 1:end÷2]
        Vn   = phase.capacity.A[2][end÷2+1:end, end÷2+1:end]
        Hₙ   = sum(diag(Vn))
        Hₙ₊₁ = sum(diag(Vn_1))

        # Compute flux through the interface using operator matrices and temperatures
        W! = phase.operator.Wꜝ[1:nx, 1:nx]
        G  = phase.operator.G[1:nx, 1:nx]
        H  = phase.operator.H[1:nx, 1:nx]
        Tₒ, Tᵧ = Tᵢ[1:nx], Tᵢ[nx+1:end]
        Interface_term = H' * W! * G * Tₒ + H' * W! * H * Tᵧ
        Interface_term = -1/(ρ*L) * sum(Interface_term)
        
        local xf_new = Hₙ + Interface_term
        err = abs(xf_new - old_xf)
        push!(residuals, err)
        println("Newton iteration $iter | xf = $xf_new | error = $err")

        # 3) Check convergence
        if err <= tol
            xf = xf_new
            break
        end
        old_xf = xf_new
        xf = xf_new

        # 4) Rebuild geometry for updated interface
        tn1, tn = tspan[2], tspan[1]
        body = Body((xx, tt, _=0) -> (xx - (xf*(tn1 - tt)/Δt + xf_new*(tt - tn)/Δt)),
                    (xx,)         -> (xx,),
                    ((0.0,1.0),),
                    false)
        # Update capacity and operator using the new body geometry
        spaceTimeMesh = phase.capacity.mesh
        spaceTimeMesh.tag = mesh.tag
        capacity = Capacity(body, spaceTimeMesh)
        operator = SpaceTimeOps(capacity.A, capacity.B, capacity.V, capacity.W, (nx, 2))
        phase = Phase(capacity, operator, (x,y,z,t)->0.0, 1.0)

        if scheme == "CN"
            s.A = build_mono_unstead_diff_moving_matrix(operator, capacity, phase.Diffusion_coeff, bc_b, bc, Δt, "CN")
            s.b = build_rhs_mono_unstead_moving_diff(operator, phase.source, capacity, bc_b, bc, Tᵢ, Δt, 0.0, "CN")
        else
            s.A = build_mono_unstead_diff_moving_matrix(operator, capacity, phase.Diffusion_coeff, bc_b, bc, Δt, "BE")
            s.b = build_rhs_mono_unstead_moving_diff(operator, phase.source, capacity, bc_b, bc, Tᵢ, Δt, 0.0, "BE")
        end
        BC_border_mono!(s.A, s.b, bc_b, capacity.mesh)
    end

    if err <= tol
        println("Newton converged after $iter iterations: xf = $xf, error = $err")
    else
        println("Newton stopped at max_iter = $max_iter: xf = $xf, error = $err")
    end

    return xf, residuals
end

function solve_MotionDiffusionUnsteadyMono1!(s::Solver, phase::Phase, Tᵢ::Vector{Float64}, Δt::Float64, Tₑ::Float64, nt::Int, bc_b::BorderConditions, bc::AbstractBoundary, body::Body, mesh::CartesianMesh, t::Vector{Float64}, scheme::String, xf; method = IterativeSolvers.gmres, kwargs...)
    if s.A === nothing
        error("Solver is not initialized. Call a solver constructor first.")
    end

    println("Solving the problem:")
    println("- Moving problem")
    println("- Monophasic problem")
    println("- Unsteady problem")
    println("- Diffusion problem")

    # Params
    nx, _ = phase.operator.size
    ρ, L = 1.0, 1.0
    max_iter, tol = 1000, 1e-15
    err = Inf
    old_xf = xf
    iter = 0
    residuals1 = Float64[]

    # First time step : Newton to compute the interface position xf1
    while (iter < max_iter) && (err > tol)
        iter += 1

        # 1) Solve the linear system
        if method == \
            A_reduced, b_reduced, rows_idx, cols_idx = remove_zero_rows_cols!(s.A, s.b)
            x_reduced = A_reduced \ b_reduced
            s.x = zeros(size(s.A, 1))
            s.x[cols_idx] .= x_reduced
        else
            log = get(kwargs, :log, false)
            if log
                s.x, s.ch = method(s.A, s.b; kwargs...)
            else
                s.x = method(s.A, s.b; kwargs...)
            end
        end
        Tᵢ = s.x

        # 2) Update volumes / compute new interface
        Vn_1 = phase.capacity.A[2][1:end÷2, 1:end÷2]
        Vn   = phase.capacity.A[2][end÷2+1:end, end÷2+1:end]
        Hₙ   = sum(diag(Vn))
        Hₙ₊₁ = sum(diag(Vn_1))

        # Compute flux
        W! = phase.operator.Wꜝ[1:nx, 1:nx]
        G = phase.operator.G[1:nx, 1:nx]
        H = phase.operator.H[1:nx, 1:nx]
        V = phase.operator.V[1:nx, 1:nx]
        Tₒ, Tᵧ = Tᵢ[1:nx], Tᵢ[nx+1:end]
        Interface_term = H' * W! * G * Tₒ + H' * W! * H * Tᵧ
        Interface_term = -1/(ρ*L) * sum(Interface_term)

        # New interface position
        res = Hₙ₊₁ - Hₙ - Interface_term
        global xf_new = Hₙ + Interface_term
        err = abs(xf_new - old_xf)
        println("Iteration $iter | xf = $xf_new | error = $err")

        # Store residuals
        push!(residuals1, err)

        # 3) Update geometry if not converged
        if err <= tol
            break
        end
        old_xf = xf_new

        # Store tn+1 and tn
        tn1 = t[2]
        tn  = t[1]

        # 4) Rebuild domain : # Add t interpolation : x - (xf*(tn1 - t)/(\Delta t) + xff*(t - tn)/(\Delta t))
        body = Body((xx,t, _=0)->(xx - (xf*(tn1 - t)/Δt + xf_new*(t - tn)/Δt)), (xx,)->(xx,), ((0.0,1.0),), false)
        spaceTimeMesh = phase.capacity.mesh
        spaceTimeMesh.tag = mesh.tag

        capacity = Capacity(body, spaceTimeMesh)
        operator = SpaceTimeOps(capacity.A, capacity.B, capacity.V, capacity.W, (nx, 2))
        phase   = Phase(capacity, operator, (x,y,z,t)->0.0, 1.0)

        if scheme == "CN"
            s.A = build_mono_unstead_diff_moving_matrix(phase.operator, phase.capacity, phase.Diffusion_coeff, bc_b, bc, Δt, "CN")
            s.b = build_rhs_mono_unstead_moving_diff(phase.operator, phase.source, phase.capacity, bc_b, bc, Tᵢ, Δt, 0.0, "CN")
        else 
            s.A = build_mono_unstead_diff_moving_matrix(phase.operator, phase.capacity, phase.Diffusion_coeff, bc_b, bc, Δt, "BE")
            s.b = build_rhs_mono_unstead_moving_diff(phase.operator, phase.source, phase.capacity, bc_b, bc, Tᵢ, Δt, 0.0, "BE")
        end
        BC_border_mono!(s.A, s.b, bc_b, phase.capacity.mesh)

    end

    if err <= tol
        println("Converged after $iter iterations with xf = $xf_new, error = $err")
    else
        println("Reached max_iter = $max_iter with xf = $xf_new, error = $err")
    end

    push!(s.states, s.x)
    println("Time : $(t[1])")
    println("Max value : $(maximum(abs.(s.x)))")

    # Reconstruct
    spaceTimeMesh = CartesianSpaceTimeMesh(mesh, t[2:3];tag=mesh.tag)
    @show xf_new
    body = Body((xx,t, _=0)->(xx - xf_new), (xx,)->(xx,), ((0.0,1.0),), false)
    capacity = Capacity(body, spaceTimeMesh)
    operator = SpaceTimeOps(capacity.A, capacity.B, capacity.V, capacity.W, (nx, 2))
    phase = Phase(capacity, operator, (x,y,z,t)->0.0, 1.0)

    s.A = build_mono_unstead_diff_moving_matrix(phase.operator, phase.capacity, phase.Diffusion_coeff, bc_b, bc, Δt, scheme)
    s.b = build_rhs_mono_unstead_moving_diff(phase.operator, phase.source, phase.capacity, bc_b, bc, s.x, Δt, t[2], scheme)
    BC_border_mono!(s.A, s.b, bc_b, capacity.mesh)

    err = Inf
    global old_xf = xf_new
    iter = 0
    residuals2 = Float64[]

    # Second time step : Newton to compute the interface position xf2
    while (iter < max_iter) && (err > tol)
        iter += 1

        # 1) Solve the linear system
        if method == \
            A_reduced, b_reduced, rows_idx, cols_idx = remove_zero_rows_cols!(s.A, s.b)
            x_reduced = A_reduced \ b_reduced
            s.x = zeros(size(s.A, 1))
            s.x[cols_idx] .= x_reduced
        else
            log = get(kwargs, :log, false)
            if log
                s.x, s.ch = method(s.A, s.b; kwargs...)
            else
                s.x = method(s.A, s.b; kwargs...)
            end
        end
        Tᵢ = s.x

        # 2) Update volumes / compute new interface
        Vn_1 = phase.capacity.A[2][1:end÷2, 1:end÷2]
        Vn   = phase.capacity.A[2][end÷2+1:end, end÷2+1:end]
        Hₙ   = sum(diag(Vn))
        Hₙ₊₁ = sum(diag(Vn_1))

        # Compute flux
        W! = phase.operator.Wꜝ[1:nx, 1:nx]
        G = phase.operator.G[1:nx, 1:nx]
        H = phase.operator.H[1:nx, 1:nx]
        V = phase.operator.V[1:nx, 1:nx]
        Tₒ, Tᵧ = Tᵢ[1:nx], Tᵢ[nx+1:end]
        Interface_term = H' * W! * G * Tₒ + H' * W! * H * Tᵧ
        Interface_term = -1/(ρ*L) * sum(Interface_term)

        # New interface position
        res = Hₙ₊₁ - Hₙ - Interface_term
        global xf_new = Hₙ + Interface_term
        err = abs(xf_new - old_xf)
        println("Iteration $iter | xf = $xf_new | error = $err")

        # Store residuals
        push!(residuals2, err)

        # 3) Update geometry if not converged
        if err <= tol
            break
        end
        old_xf = xf_new

        # Store tn+1 and tn
        tn1 = t[3]
        tn  = t[2]

        # 4) Rebuild domain : # Add t interpolation : x - (xf*(tn1 - t)/(\Delta t) + xff*(t - tn)/(\Delta t))
        body = Body((xx,t, _=0)->(xx - (xf*(tn1 - t)/Δt + xf_new*(t - tn)/Δt)), (xx,)->(xx,), ((0.0,1.0),), false)
        spaceTimeMesh = phase.capacity.mesh
        spaceTimeMesh.tag = mesh.tag
        
        capacity = Capacity(body, spaceTimeMesh)
        operator = SpaceTimeOps(capacity.A, capacity.B, capacity.V, capacity.W, (nx, 2))
        phase   = Phase(capacity, operator, (x,y,z,t)->0.0, 1.0)

        if scheme == "CN"
            s.A = build_mono_unstead_diff_moving_matrix(phase.operator, phase.capacity, phase.Diffusion_coeff, bc_b, bc, Δt, "CN")
            s.b = build_rhs_mono_unstead_moving_diff(phase.operator, phase.source, phase.capacity, bc_b, bc, Tᵢ, Δt, 0.0, "CN")
        else 
            s.A = build_mono_unstead_diff_moving_matrix(phase.operator, phase.capacity, phase.Diffusion_coeff, bc_b, bc, Δt, "BE")
            s.b = build_rhs_mono_unstead_moving_diff(phase.operator, phase.source, phase.capacity, bc_b, bc, Tᵢ, Δt, 0.0, "BE")
        end

        BC_border_mono!(s.A, s.b, bc_b, phase.capacity.mesh)
    end

    if err <= tol
        println("Converged after $iter iterations with xf = $xf_new, error = $err")
    else
        println("Reached max_iter = $max_iter with xf = $xf_new, error = $err")
    end

    push!(s.states, s.x)
    println("Time : $(t[2])")
    println("Max value : $(maximum(abs.(s.x)))")

    return residuals1, residuals2
end

residuals1, residuals2 = solve_MotionDiffusionUnsteadyMono1!(s, Fluide, u0, Δt, Tend, nt, bc_b, bc, body, mesh, t, "BE", xf)

# Plot the residuals
fig = Figure()
ax = Axis(fig[1, 1], title = "Residuals", xlabel = "Iteration", ylabel = "Residuals")
lines!(ax, 1:length(residuals1), log10.(residuals1), color = :blue, label = "Residuals 1")
lines!(ax, 1:length(residuals2), log10.(residuals2), color = :red, label = "Residuals 2")
axislegend(ax, position = :rt)
display(fig)

