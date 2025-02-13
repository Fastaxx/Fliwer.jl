using Fliwer
using IterativeSolvers
using LinearAlgebra
using SparseArrays
using SpecialFunctions, LsqFit
using CairoMakie

### 1D Test Case : Monophasic Unsteady Diffusion Equation inside a moving body
# Define the spatial mesh
nx = 80
lx = 10.
x0 = 0.
domain = ((x0, lx),)
mesh = CartesianMesh((nx,), (lx,), (x0,))
x = range(x0, stop = lx, length = nx+1)

Δt = 0.25 * (lx/nx)^2
Tend = 1.0
nt = 10
t = [i*Δt for i in 0:nt]

# Define the body
xf = 0.03*lx   # Interface position$
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
function stefan_1d_1ph_analytical(x::Float64)
    t = Tend
    λ = 0.029725478966757847
    return 1.0 - 1.0/erf(λ) * erf(x/(2*sqrt(t)))
end
x = range(x0, stop=lx, length=nx+1)
u0ₒ = [stefan_1d_1ph_analytical(x[i]-xf) for i in 1:nx+1]
u0ᵧ = zeros((nx+1))
u0 = vcat(u0ₒ, u0ᵧ)

# Define the solver
s = MotionDiffusionUnsteadyMono(Fluide, bc_b, bc, Δt, Tend, u0, "BE")

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
    max_iter, tol = 1000, 1e-10
    err = Inf
    old_xf = xf
    iter = 0
    residuals1 = Float64[]
    xf_log = Float64[]

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
        @show Interface_term
        Interface_term = -1/(ρ*L) * sum(Interface_term)

        # New interface position
        res = Hₙ₊₁ - Hₙ - Interface_term
        α = 0.5
        global xf_new = xf + α*res
        err = abs(xf_new - old_xf)
        println("Iteration $iter | xf = $xf_new | error = $err | res = $res")

        # Store residuals
        push!(residuals1, err)

        # 3) Update geometry if not converged
        if err <= tol
            push!(xf_log, xf_new)
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
    body = Body((xx,t, _=0)->(xx - xf_new), (xx,)->(xx,), ((0.0,1.0),), false) # Possibly add previous velocity to have a sdf with t.
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
        global xf_new = old_xf + res
        err = abs(xf_new - old_xf)
        println("Iteration $iter | xf = $xf_new | error = $err | res = $res")

        # Store residuals
        push!(residuals2, err)

        # 3) Update geometry if not converged
        if err <= tol
            push!(xf_log, xf_new)
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

    # Reconstruct
    spaceTimeMesh = CartesianSpaceTimeMesh(mesh, t[3:4];tag=mesh.tag)
    body = Body((xx,t, _=0)->(xx - xf_new), (xx,)->(xx,), ((0.0,1.0),), false) # Possibly add previous velocity to have a sdf with t.
    capacity = Capacity(body, spaceTimeMesh)
    operator = SpaceTimeOps(capacity.A, capacity.B, capacity.V, capacity.W, (nx, 2))
    phase = Phase(capacity, operator, (x,y,z,t)->0.0, 1.0)

    s.A = build_mono_unstead_diff_moving_matrix(phase.operator, phase.capacity, phase.Diffusion_coeff, bc_b, bc, Δt, scheme)
    s.b = build_rhs_mono_unstead_moving_diff(phase.operator, phase.source, phase.capacity, bc_b, bc, s.x, Δt, t[3], scheme)
    BC_border_mono!(s.A, s.b, bc_b, capacity.mesh)

    err = Inf
    global old_xf = xf_new
    iter = 0
    residuals3 = Float64[]

    # Third time step : Newton to compute the interface position xf3
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
        global xf_new = old_xf + res
        err = abs(xf_new - old_xf)
        println("Iteration $iter | xf = $xf_new | error = $err | res = $res")

        # Store residuals
        push!(residuals3, err)

        # 3) Update geometry if not converged
        if err <= tol
            push!(xf_log, xf_new)
            break
        end
        old_xf = xf_new

        # Store tn+1 and tn
        tn1 = t[4]
        tn  = t[3]

        # 4) Rebuild domain : # Add t interpolation : x - (xf*(tn1 - t)/(\Delta t) + xff*(t - tn)/(\Delta t))
        body = Body((xx,t, _=0)->(xx - (xf_log[end]*(tn1 - t)/Δt + xf_new*(t - tn)/Δt)), (xx,)->(xx,), ((0.0,1.0),), false)
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
    println("Time : $(t[3])")
    println("Max value : $(maximum(abs.(s.x)))")

    # Reconstruct
    spaceTimeMesh = CartesianSpaceTimeMesh(mesh, t[4:5];tag=mesh.tag)
    body = Body((xx,t, _=0)->(xx - xf_new), (xx,)->(xx,), ((0.0,1.0),), false) # Possibly add previous velocity to have a sdf with t.
    capacity = Capacity(body, spaceTimeMesh)
    operator = SpaceTimeOps(capacity.A, capacity.B, capacity.V, capacity.W, (nx, 2))
    phase = Phase(capacity, operator, (x,y,z,t)->0.0, 1.0)

    s.A = build_mono_unstead_diff_moving_matrix(phase.operator, phase.capacity, phase.Diffusion_coeff, bc_b, bc, Δt, scheme)
    s.b = build_rhs_mono_unstead_moving_diff(phase.operator, phase.source, phase.capacity, bc_b, bc, s.x, Δt, t[4], scheme)

    BC_border_mono!(s.A, s.b, bc_b, capacity.mesh)

    err = Inf
    global old_xf = xf_new
    iter = 0
    residuals4 = Float64[]

    # Fourth time step : Newton to compute the interface position xf4
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
        global xf_new = old_xf + res
        err = abs(xf_new - old_xf)
        println("Iteration $iter | xf = $xf_new | error = $err | res = $res")

        # Store residuals
        push!(residuals4, err)

        # 3) Update geometry if not converged
        if err <= tol
            push!(xf_log, xf_new)
            break
        end
        old_xf = xf_new

        # Store tn+1 and tn
        tn1 = t[5]
        tn  = t[4]

        # 4) Rebuild domain : # Add t interpolation : x - (xf*(tn1 - t)/(\Delta t) + xff*(t - tn)/(\Delta t))
        body = Body((xx,t, _=0)->(xx - (xf_log[end]*(tn1 - t)/Δt + xf_new*(t - tn)/Δt)), (xx,)->(xx,), ((0.0,1.0),), false)
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
    println("Time : $(t[4])")
    println("Max value : $(maximum(abs.(s.x)))")


    
    return residuals1, residuals2, residuals3, residuals4, xf_log
end

function solve_MotionDiffusionUnsteadyMonoN!(s::Solver, phase::Phase, Tᵢ::Vector{Float64},
    Δt::Float64, Tₑ::Float64, nt::Int, bc_b::BorderConditions, bc::AbstractBoundary,
    body::Body, mesh::CartesianMesh, t::Vector{Float64}, scheme::String, xf; 
    method = IterativeSolvers.gmres, kwargs...)

    # Check solver initialization
    if s.A === nothing
        error("Solver is not initialized. Call a solver constructor first.")
    end

    println("Solving the problem:")
    println("- Moving problem")
    println("- Monophasic problem")
    println("- Unsteady problem")
    println("- Diffusion problem")

    nx, _ = phase.operator.size
    ρ, L = 1.0, 1.0
    max_iter = 1000
    tol      = 1e-10

    # Log residuals and interface positions for each time step:
    residuals = [Float64[] for i in 1:nt]
    xf_log = Float64[]

    # Initialize the interface position for the first time step:
    current_xf = xf

    # Loop over each time step (using t[i] and t[i+1])
    for step in 1:nt
        println("\n===== Time Step $step: t = $(t[step]) → $(t[step+1]) =====")
        err = Inf
        iter = 0

        # Initialize new_xf in the outer scope of the Newton loop
        new_xf = current_xf
        xf = current_xf

        # Newton iterations for this time step:
        while (iter < max_iter) && (err > tol)
            iter += 1

            # --- 1) Solve the linear system ---
            # (The block below reproduces your two alternatives.)
            if method == :reduced  # if you want to use the reduced system
                A_reduced, b_reduced, rows_idx, cols_idx = remove_zero_rows_cols!(s.A, s.b)
                x_reduced = A_reduced \ b_reduced
                s.x = zeros(size(s.A,1))
                s.x[cols_idx] .= x_reduced
            else
                log_flag = get(kwargs, :log, false)
                if log_flag
                    s.x, s.ch = method(s.A, s.b; kwargs...)
                else
                    s.x = method(s.A, s.b; kwargs...)
                end
            end
            # Update the current solution
            Tᵢ = s.x

            # --- 2) Compute the interface update ---
            # Here we “extract” the two sub-blocks of the capacity matrix
            Vn_1 = phase.capacity.A[2][1:end÷2, 1:end÷2]
            Vn   = phase.capacity.A[2][end÷2+1:end, end÷2+1:end]
            Hₙ   = sum(diag(Vn))
            Hₙ₊₁ = sum(diag(Vn_1))

            # Compute the flux (using the operators on the first nx rows/cols)
            W_part = phase.operator.Wꜝ[1:nx, 1:nx]
            G      = phase.operator.G[1:nx, 1:nx]
            H      = phase.operator.H[1:nx, 1:nx]
            Tₒ, Tᵧ = s.x[1:nx], s.x[nx+1:end]
            Interface_term = H' * W_part * G * Tₒ + H' * W_part * H * Tᵧ
            Interface_term = -1/(ρ*L) * sum(Interface_term)

            # The residual (difference between new and old volumes minus the flux term)
            res_val = Hₙ₊₁ - Hₙ - Interface_term

            # Update the interface position:
            if step == 1
                α = 0.5
                new_xf = current_xf + α*res_val
            else
                new_xf = current_xf + res_val
            end

            err = abs(new_xf - current_xf)
            println("  Newton iter $iter: xf = $new_xf, error = $err, res = $res_val")
            push!(residuals[step], err)

            # If converged for this time step, record the new interface position and exit Newton loop.
            if err <= tol
                push!(xf_log, new_xf)
                break
            end

            # --- 3) Update geometry (if not converged) ---
            # Set the time points for interpolation on this time interval:
            tn  = t[step]
            tn1 = t[step+1]
            # The interpolation function here uses the old and new interface positions:
            body = Body((xx, t_val, _=0) -> (xx - (xf*(tn1 - t_val)/Δt + new_xf*(t_val - tn)/Δt)),
                        (xx,) -> xx,
                        ((0.0, 1.0),), false)
            # Rebuild the space-time mesh (keeping the same tag as the original mesh)
            spaceTimeMesh = phase.capacity.mesh
            spaceTimeMesh.tag = mesh.tag

            # Rebuild the capacity, operator, and phase:
            capacity = Capacity(body, spaceTimeMesh)
            operator = SpaceTimeOps(capacity.A, capacity.B, capacity.V, capacity.W, (nx, 2))
            phase = Phase(capacity, operator, (x,y,z,t_val)->0.0, 1.0)

            # --- 4) Rebuild the linear system for the next Newton iteration ---
            if scheme == "CN"
                s.A = build_mono_unstead_diff_moving_matrix(phase.operator, phase.capacity,
                        phase.Diffusion_coeff, bc_b, bc, Δt, "CN")
                s.b = build_rhs_mono_unstead_moving_diff(phase.operator, phase.source, phase.capacity,
                        bc_b, bc, s.x, Δt, 0.0, "CN")
            else
                s.A = build_mono_unstead_diff_moving_matrix(phase.operator, phase.capacity,
                        phase.Diffusion_coeff, bc_b, bc, Δt, "BE")
                s.b = build_rhs_mono_unstead_moving_diff(phase.operator, phase.source, phase.capacity,
                        bc_b, bc, s.x, Δt, 0.0, "BE")
            end
            BC_border_mono!(s.A, s.b, bc_b, phase.capacity.mesh)

            # Update the interface for the next Newton iteration
            current_xf = new_xf
        end

        if err <= tol
            println("Time step $step converged after $iter iterations with xf = $new_xf, error = $err")
        else
            println("Time step $step reached max_iter = $max_iter with xf = $new_xf, error = $err")
        end

        # --- 5) Save the state and report information ---
        push!(s.states, s.x)
        println("  Time: $(t[step]), Max value in s.x: $(maximum(abs.(s.x)))")

        # --- 6) Reconstruct the domain for the next time step (if any) ---
        if step < nt
            # Here we build a new space-time mesh using the next time interval [t[step+1], t[step+2]]
            # (Note: for the last step there is no “next” time interval.)
            spaceTimeMesh = CartesianSpaceTimeMesh(mesh, t[step+1:step+2]; tag=mesh.tag)
            body = Body((xx, t_val, _=0) -> (xx - new_xf),
                        (xx,) -> xx,
                        ((0.0, 1.0),), false)
            capacity = Capacity(body, spaceTimeMesh)
            operator = SpaceTimeOps(capacity.A, capacity.B, capacity.V, capacity.W, (nx, 2))
            phase = Phase(capacity, operator, (x,y,z,t_val)->0.0, 1.0)

            # Update s.A and s.b for the next time step
            s.A = build_mono_unstead_diff_moving_matrix(phase.operator, phase.capacity,
                        phase.Diffusion_coeff, bc_b, bc, Δt, scheme)
            s.b = build_rhs_mono_unstead_moving_diff(phase.operator, phase.source, phase.capacity,
                        bc_b, bc, s.x, Δt, t[step+1], scheme)
            BC_border_mono!(s.A, s.b, bc_b, capacity.mesh)
        end
    end

    return residuals, xf_log
end

# Solve the problem
residuals, xf_log = solve_MotionDiffusionUnsteadyMonoN!(s, Fluide, u0, Δt, Tend, nt, bc_b, bc, body, mesh, t, "BE", xf)

#residuals1, residuals2, residuals3, residuals4, xf_log = solve_MotionDiffusionUnsteadyMono1!(s, Fluide, u0, Δt, Tend, nt, bc_b, bc, body, mesh, t, "BE", xf)

# Plot the residuals
fig = Figure()
ax = Axis(fig[1, 1], title = "Residuals", xlabel = "Iteration", ylabel = "Residuals")
for i in 1:nt
    lines!(ax, 1:length(residuals[i]), log10.(residuals[i]), label = "Residuals $i")
end
axislegend(ax, position = :rt)
display(fig)
readline()
# Plot the interface position
fig = Figure()
ax = Axis(fig[1, 1], title = "Interface position", xlabel = "Iteration", ylabel = "Interface position")
lines!(ax, 1:length(xf_log), xf_log, color = :blue, label = "Interface position")
axislegend(ax, position = :rt)
display(fig)

"""
# Plot the residuals
fig = Figure()
ax = Axis(fig[1, 1], title = "Residuals", xlabel = "Iteration", ylabel = "Residuals")
lines!(ax, 1:length(residuals1), log10.(residuals1), color = :blue, label = "Residuals 1")
lines!(ax, 1:length(residuals2), log10.(residuals2), color = :red, label = "Residuals 2")
lines!(ax, 1:length(residuals3), log10.(residuals3), color = :green, label = "Residuals 3")
lines!(ax, 1:length(residuals4), log10.(residuals4), color = :black, label = "Residuals 4")
axislegend(ax, position = :rt)
display(fig)
readline()

# Plot the interface position
fig = Figure()
ax = Axis(fig[1, 1], title = "Interface position", xlabel = "Iteration", ylabel = "Interface position")
lines!(ax, 1:length(xf_log), xf_log, color = :blue, label = "Interface position")
axislegend(ax, position = :rt)
display(fig)
readline()

# Plot the solution
plot_solution(s, mesh, body, capacity; state_i = 1)
"""