function DarcyFlow(phase::Phase, bc_b::BorderConditions, bc_i::AbstractBoundary)
    println("Création du solveur:")
    println("- Darcy Flow")
    println("- Steady problem")
    println("- Monophasic")

    s = Solver(Steady, Monophasic, Diffusion, nothing, nothing, nothing, ConvergenceHistory(), [])

    s.A = build_mono_stead_diff_matrix(phase.operator, phase.capacity, phase.Diffusion_coeff, bc_b, bc_i)
    s.b = build_rhs_mono_stead_diff(phase.operator, phase.source, phase.capacity, bc_b, bc_i)

    BC_border_mono!(s.A, s.b, bc_b, phase.capacity.mesh)

    return s
end

function solve_DarcyFlow!(s::Solver, phase::Phase; method::Function = gmres, kwargs...)
    if s.A === nothing
        error("Solver is not initialized. Call a solver constructor first.")
    end

    n = Int(size(s.A, 1) / 2)

    if method == \
        # Remove zero rows and columns for direct solver
        s.A, s.b, rows_idx, cols_idx = remove_zero_rows_cols!(s.A, s.b)
        # Solve the reduced system
        x_reduced = s.A \ s.b
        # Reconstruct the full solution vector
        s.x = zeros(2n)
        s.x[cols_idx] = x_reduced
    else
        # Use iterative solver directly
        kwargs_nt = (; kwargs...)
        log = get(kwargs_nt, :log, false)
        if log
            s.x, s.ch = method(s.A, s.b; kwargs...)
        else
            s.x = method(s.A, s.b; kwargs...)
        end
    end

    push!(s.states, s.x)
end

function solve_darcy_velocity(solver, Fluide; state_i=1)
    cell_types = Fluide.capacity.cell_types
    pₒ = solver.states[state_i][1:div(end,2)]
    pᵧ = solver.states[state_i][div(end,2)+1:end]

    pₒ[cell_types .== 0] .= NaN
    pᵧ[cell_types .== 0] .= NaN
    pᵧ[cell_types .== 1] .= NaN

    p = vcat(pₒ, pᵧ)

    # Compute the velocity field
    u = - ∇(Fluide.operator, p)
    return u
end


# Unsteady Darcy Flow

function DarcyFlowUnsteady(phase::Phase, bc_b::BorderConditions, bc_i::AbstractBoundary, Δt::Float64, Tend::Float64, Tᵢ::Vector{Float64})
    println("Création du solveur:")
    println("- Darcy Flow")
    println("- Unsteady problem")
    println("- Monophasic")

    s = Solver(Unsteady, Monophasic, Diffusion, nothing, nothing, nothing, ConvergenceHistory(), [])

    s.A = build_mono_unstead_diff_matrix(phase.operator, phase.capacity, phase.Diffusion_coeff, bc_b, bc_i, Δt)
    s.b = build_rhs_mono_unstead_diff(phase.operator, phase.source, phase.capacity, bc_b, bc_i, Tᵢ, Δt, 0.0)
    return s
end

function solve_DarcyFlowUnsteady!(s::Solver, phase::Phase, Tᵢ::Vector{Float64}, Δt::Float64, Tₑ::Float64, bc_b::BorderConditions, bc_i::AbstractBoundary; method::Function = gmres, kwargs...)
    if s.A === nothing
        error("Solver is not initialized. Call a solver constructor first.")
    end

    n = Int(size(s.A, 1) / 2)

    t=0.0
    while t<Tₑ
        t+=Δt
        println("Time: ", t)
        s.b = build_rhs_mono_unstead_diff(phase.operator, phase.source, phase.capacity, bc_b, bc_i, Tᵢ, Δt, t)
        BC_border_mono!(s.A, s.b, bc_b, phase.capacity.mesh)

        if method == \
            # Remove zero rows and columns for direct solver
            A_reduced, b_reduced, _, cols_idx = remove_zero_rows_cols!(s.A, s.b)
            # Solve the reduced system
            x_reduced = A_reduced \ b_reduced
            # Reconstruct the full solution vector
            s.x = zeros(2n)
            s.x[cols_idx] = x_reduced
        else
            kwargs_nt = (; kwargs...)
            log = get(kwargs_nt, :log, false)
            if log
                s.x, s.ch = method(s.A, s.b; kwargs...)
            else
                s.x = method(s.A, s.b; kwargs...)
            end
        end

        push!(s.states, s.x)
        @show maximum(s.x)
        Tᵢ = s.x
    end
end