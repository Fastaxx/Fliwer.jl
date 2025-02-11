function MotionDiffusionUnsteadyMono(phase::Phase, bc_b::BorderConditions, bc_i::AbstractBoundary, Δt::Float64, Tₑ::Float64, Tᵢ::Vector{Float64}, scheme::String)
    println("Création du solveur:")
    println("- Moving problem")
    println("- Monophasic problem")
    println("- Unsteady problem")
    println("- Diffusion problem")
    
    s = Solver(Unsteady, Monophasic, Diffusion, nothing, nothing, nothing, ConvergenceHistory(), [])
    
    if scheme == "CN"
        s.A = build_mono_unstead_diff_moving_matrix(phase.operator, phase.capacity, phase.Diffusion_coeff, bc_b, bc_i, Δt, "CN")
        s.b = build_rhs_mono_unstead_moving_diff(phase.operator, phase.source, phase.capacity, bc_b, bc_i, Tᵢ, Δt, 0.0, "CN")
    else 
        s.A = build_mono_unstead_diff_moving_matrix(phase.operator, phase.capacity, phase.Diffusion_coeff, bc_b, bc_i, Δt, "BE")
        s.b = build_rhs_mono_unstead_moving_diff(phase.operator, phase.source, phase.capacity, bc_b, bc_i, Tᵢ, Δt, 0.0, "BE")
    end
    BC_border_mono!(s.A, s.b, bc_b, phase.capacity.mesh)

    return s
end

function solve_MotionDiffusionUnsteadyMono!(s::Solver, phase::Phase, Tᵢ::Vector{Float64}, Δt::Float64, Tₑ::Float64, nt::Int, bc_b::BorderConditions, bc::AbstractBoundary, body::Body, mesh::CartesianMesh, t::Vector{Float64}, scheme::String, xf; method = IterativeSolvers.gmres, kwargs...)
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
    max_iter, tol = 100, 1e-20
    err = Inf
    old_xf = xf
    iter = 0
    residuals = Float64[]

    # First time step : Newton 
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
        xf_new = Hₙ + Interface_term
        err = abs(xf_new - old_xf)
        println("Iteration $iter | xf = $xf_new | error = $err")

        # Store residuals
        push!(residuals, err)

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
        Fluide   = Phase(capacity, operator, (x,y,z,t)->0.0, 1.0)

        s = MotionDiffusionUnsteadyMono(Fluide, bc_b, bc, Δt, Tₑ, Tᵢ, scheme)

    end

    if err <= tol
        println("Converged after $iter iterations with xf = $xf, error = $err")
    else
        println("Reached max_iter = $max_iter with xf = $xf, error = $err")
    end

    push!(s.states, s.x)
    println("Time : $(t[1])")
    println("Max value : $(maximum(abs.(s.x)))")
    
    for i in 2 :nt
        println("Time : $(t[i])")
        spaceTimeMesh = CartesianSpaceTimeMesh(mesh, t[i:i+1];tag=mesh.tag)
        capacity = Capacity(body, spaceTimeMesh)
        operator = SpaceTimeOps(capacity.A, capacity.B, capacity.V, capacity.W, (nx, 2))

        s.A = build_mono_unstead_diff_moving_matrix(operator, capacity, phase.Diffusion_coeff, bc_b, bc, Δt, scheme)
        s.b = build_rhs_mono_unstead_moving_diff(operator, phase.source, capacity, bc_b, bc, Tᵢ, Δt, t[i], scheme)
        BC_border_mono!(s.A, s.b, bc_b, capacity.mesh)

        # CFL log
        Vn_1 = capacity.A[2][1:end÷2, 1:end÷2]
        Vn = capacity.A[2][end÷2+1:end, end÷2+1:end]

        # Solve system
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

        push!(s.states, s.x)
        @show maximum(abs.(s.x))
        Tᵢ = s.x
    end

end


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
    residuals = Float64[]  # Store residual at each iteration
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
        Interface_term = -1/(ρ*L) * sum(Interface_term)
        @show Interface_term
        # New interface position
        res = Hₙ₊₁ - Hₙ - Interface_term
        @show -diag(Vn_1 - Vn) * ρ * L - H' * W! * G * Tₒ + H' * W! * H * Tᵧ
        @show res
        @show Hₙ₊₁ - Hₙ
        xf_new = Hₙ + Interface_term
        
        err = abs(xf_new - old_xf)
        println("Iteration $iter | xf = $xf_new | error = $err")

        # Store residuals
        push!(residuals, err)

        # 3) Update geometry if not converged
        if err <= tol
            break
        end
        old_xf = xf_new

        # Store tn+1 and tn
        tn1 = t[2]
        tn  = t[1]

        # 4) Rebuild domain : # Add t interpolation : x - (xf*(tn1 - t)/(\Delta t) + xff*(t - tn)/(\Delta t))
        body = Body((xx,t, _=0)->(xx - (xf*(tn1 - t)/Δt + xf_new*(t - tn)/Δt)), (xx,)->(xx,), domain, false)
        #identify!(mesh, body)
        spaceTimeMesh.tag = mesh.tag

        capacity = Capacity(body, spaceTimeMesh)
        operator = SpaceTimeOps(capacity.A, capacity.B, capacity.V, capacity.W, (nx, 2))
        Fluide   = Phase(capacity, operator, (x,y,z,t)->0.0, 1.0)

        s = MotionDiffusionUnsteadyMono(Fluide, bc_b, bc, Δt, Tend, u0, "BE")

        # recompute the linear system
        A_reduced, b_reduced, rows_idx, cols_idx = remove_zero_rows_cols!(s.A, s.b)
        x_reduced = A_reduced \ b_reduced
        s.x = zeros(size(s.A, 1))
        s.x[cols_idx] .= x_reduced
        Tᵢ = s.x
    end

    if err <= tol
        println("Converged after $iter iterations with xf = $xf, error = $err")
    else
        println("Reached max_iter = $max_iter with xf = $xf, error = $err")
    end

    return xf,residuals
end

