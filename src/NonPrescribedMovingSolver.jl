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

    return s
end

function solve_MotionDiffusionUnsteadyMono!(s::Solver, phase::Phase, Tᵢ::Vector{Float64}, Δt::Float64, Tₑ::Float64, nt::Int, bc_b::BorderConditions, bc::AbstractBoundary, body::Body, mesh::CartesianMesh, t::Vector{Float64}, scheme::String; method = IterativeSolvers.gmres, kwargs...)
    if s.A === nothing
        error("Solver is not initialized. Call a solver constructor first.")
    end

    println("Solving the problem:")
    println("- Moving problem")
    println("- Monophasic problem")
    println("- Unsteady problem")
    println("- Diffusion problem")

    nx, _ = phase.operator.size
    # Solve system for the initial condition
    if method == \
        A_reduced, b_reduced, rows_idx, cols_idx = remove_zero_rows_cols!(s.A, s.b)
        # Compute condition number
        #cnum = cond(Array(A_reduced), 2)
        cnum = 0.0
        push!(cond_log, cnum)
        push!(minV_log, minimum(x for x in phase.capacity.V if x != 0))
        push!(maxV_log, maximum(phase.capacity.V))
        push!(minW_log, minimum(x for x in phase.capacity.W[1] if x != 0))
        push!(maxW_log, maximum(phase.capacity.W[1]))
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

    # Newton initial
    Vn_1 = phase.capacity.A[2][1:end÷2, 1:end÷2]
    Vn = phase.capacity.A[2][end÷2+1:end, end÷2+1:end]

    # Height function Hn and Hn+1,k : Sum of the Volume of the cells
    Hⁿ = sum(diag(Vn))
    Hⁿ⁺¹ = sum(diag(Vn_1))
    ΔH = Hⁿ⁺¹ - Hⁿ
    @show Hⁿ, Hⁿ⁺¹, ΔH

    # Interface term : 
    W! = phase.operator.Wꜝ[1:nx, 1:nx]
    G = phase.operator.G[1:nx, 1:nx]
    H = phase.operator.H[1:nx, 1:nx]
    V = phase.operator.V[1:nx, 1:nx]
    Tₒ, Tᵧ = Tᵢ[1:nx], Tᵢ[nx+1:end]
    Interface_term = H' * W! * G * Tₒ + H' * W! * H * Tᵧ
    Interface_term = -sum(Interface_term)
    @show Interface_term
    
    # Update the height function
    res = Hⁿ⁺¹ - (Hⁿ + Interface_term)
    Hⁿ⁺¹ = Hⁿ + Interface_term

    @show res
    @show Hⁿ⁺¹
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