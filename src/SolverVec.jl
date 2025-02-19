"""
    mutable struct SolverVec{TT<:TimeType, PT<:PhaseType, ET<:EquationType}

The `SolverVec` struct represents a solver for a specific type of problem - Vector problem.

# Fields
- `time_type::TT`: The type of time used in the solver : `Steady` or `Unsteady`.
- `phase_type::PT`: The type of phase used in the solver : `Monophasic` or `Diphasic`.
- `equation_type::ET`: The type of equation used in the solver : `Diffusion`, `Advection` or `DiffusionAdvection`.
- `A::Union{SparseMatrixCSC{Float64, Int}, Nothing}`: The coefficient matrix A of the equation system, if applicable.
- `b::Union{Vector{Float64}, Nothing}`: The right-hand side vector b of the equation system, if applicable.
- `x::Union{Vector{Float64}, Nothing}`: The solution vector x of the equation system, if applicable.
- `states::Vector{Any}`: The states of the system at different times, if applicable.

"""
mutable struct SolverVec{TT<:TimeType, PT<:PhaseType, ET<:EquationType}
    time_type::TT
    phase_type::PT
    equation_type::ET
    A::Union{SparseMatrixCSC{Float64, Int}, Nothing}
    b::Union{Vector{Float64}, Nothing}
    x::Union{Vector{Float64}, Nothing}
    ch::IterativeSolvers.ConvergenceHistory
    states::Vector{Any}
end

# Vector Diffusion - Unsteady - Monophasic
function DiffusionVecUnsteadyMono(phase::VectorPhase{1}, bc::Tuple{BorderConditions}, ic::Tuple{AbstractBoundary}, Δt::Float64, Tend::Float64, u0)
    println("Création du solveur:")
    println("- Vector problem")
    println("- Monophasic problem")
    println("- Unsteady problem")
    println("- Diffusion problem")
    
    s = SolverVec(Unsteady, Monophasic, Diffusion, nothing, nothing, nothing, ConvergenceHistory(), [])

    A_u = build_mono_unstead_diff_vec_matrix(phase, bc, ic, Δt)
    b_u = build_mono_unstead_diff_vec_rhs(phase, bc, ic, Δt, u0, 0.0)

    s.A = A_u
    s.b = b_u
    # BC_border_mono!(s, phase, bc)

    return s
end

function DiffusionVecUnsteadyMono(phase::VectorPhase{2}, bc::Tuple{BorderConditions, BorderConditions}, ic::Tuple{AbstractBoundary, AbstractBoundary}, Δt::Float64, Tend::Float64, u0x, u0y)
    println("Création du solveur:")
    println("- Vector problem")
    println("- Monophasic problem")
    println("- Unsteady problem")
    println("- Diffusion problem")
    
    s = SolverVec(Unsteady, Monophasic, Diffusion, nothing, nothing, nothing, ConvergenceHistory(), [])

    A_u, A_v = build_mono_unstead_diff_vec_matrix(phase, bc, ic, Δt)
    b_u, b_v = build_mono_unstead_diff_vec_rhs(phase, bc, ic, Δt, u0x, u0y, 0.0)

    s.A = blockdiag(A_u, A_v)
    s.b = vcat(b_u, b_v)
    # BC_border_mono!(s, phase, bc)

    return s
end

function DiffusionVecUnsteadyMono(phase::VectorPhase{3}, bc::Tuple{BorderConditions, BorderConditions, BorderConditions}, ic::Tuple{AbstractBoundary, AbstractBoundary, AbstractBoundary}, Δt::Float64, Tend::Float64, u0x, u0y, u0z)
    println("Création du solveur:")
    println("- Vector problem")
    println("- Monophasic problem")
    println("- Unsteady problem")
    println("- Diffusion problem")
    
    s = SolverVec(Unsteady, Monophasic, Diffusion, nothing, nothing, nothing, ConvergenceHistory(), [])

    A_u, A_v, A_w = build_mono_unstead_diff_vec_matrix(phase, bc, ic, Δt)
    b_u, b_v, b_w = build_mono_unstead_diff_vec_rhs(phase, bc, ic, Δt, u0x, u0y, u0z, 0.0)

    s.A = blockdiag(A_u, A_v, A_w)
    s.b = vcat(b_u, b_v, b_w)
    # BC_border_mono!(s, phase, bc)

    return s
end

function build_mono_unstead_diff_vec_matrix(phase::VectorPhase{1}, bc, ic, Δt)
    operator = phase.operator[1]
    capacity = phase.capacity[1]
    μ = phase.Diffusion_coeff

    # For u-component
    A_u = build_mono_unstead_diff_matrix(operator, capacity, μ, bc[1], ic[1], Δt)

    return A_u
end

function build_mono_unstead_diff_vec_matrix(phase::VectorPhase{2}, bc, ic, Δt)
    operator_u, operator_v = phase.operator
    capacity_u, capacity_v = phase.capacity
    μ = phase.Diffusion_coeff

    # For u-component
    A_u = build_mono_unstead_diff_matrix(operator_u, capacity_u, μ, bc[1], ic[1], Δt)

    # For v-component
    A_v = build_mono_unstead_diff_matrix(operator_v, capacity_v, μ, bc[2], ic[2], Δt)

    return A_u, A_v
end

function build_mono_unstead_diff_vec_matrix(phase::VectorPhase{3}, bc, ic, Δt)
    operator_u, operator_v, operator_w = phase.operator
    capacity_u, capacity_v, capacity_w = phase.capacity
    μ = phase.Diffusion_coeff

    # For u-component
    A_u = build_mono_unstead_diff_matrix(operator_u, capacity_u, μ, bc[1], ic[1], Δt)

    # For v-component
    A_v = build_mono_unstead_diff_matrix(operator_v, capacity_v, μ, bc[2], ic[2], Δt)

    # For w-component
    A_w = build_mono_unstead_diff_matrix(operator_w, capacity_w, μ, bc[3], ic[3], Δt)

    return A_u, A_v, A_w
end

function build_mono_unstead_diff_vec_rhs(phase::VectorPhase{1}, bc, ic, Δt, u0, t)
    operator = phase.operator[1]
    capacity = phase.capacity[1]
    f = phase.source[1]
    μ = phase.Diffusion_coeff

    # For u-component
    b_u = build_rhs_mono_unstead_diff(operator, f, capacity, bc[1], ic[1], vcat(u0[1], u0[2]), Δt, t)

    return b_u
end

function build_mono_unstead_diff_vec_rhs(phase::VectorPhase{2}, bc, ic, Δt, u0x, u0y, t)
    operator_u, operator_v = phase.operator
    capacity_u, capacity_v = phase.capacity
    fu, fv = phase.source
    μ = phase.Diffusion_coeff

    # For u-component
    b_u = build_rhs_mono_unstead_diff(operator_u, fu, capacity_u, bc[1], ic[1], vcat(u0x[1], u0x[2]), Δt, t)

    # For v-component
    b_v = build_rhs_mono_unstead_diff(operator_v, fv, capacity_v, bc[2], ic[2], vcat(u0y[1], u0y[2]), Δt, t)

    return b_u, b_v
end

function build_mono_unstead_diff_vec_rhs(phase::VectorPhase{3}, bc, ic, Δt, u0x, u0y, u0z, t)
    operator_u, operator_v, operator_w = phase.operator
    capacity_u, capacity_v, capacity_w = phase.capacity
    fu, fv, fw = phase.source
    μ = phase.Diffusion_coeff

    # For u-component
    b_u = build_rhs_mono_unstead_diff(operator_u, fu, capacity_u, bc[1], ic[1], vcat(u0x[1], u0x[2]), Δt, t)

    # For v-component
    b_v = build_rhs_mono_unstead_diff(operator_v, fv, capacity_v, bc[2], ic[2], vcat(u0y[1], u0y[2]), Δt, t)

    # For w-component
    b_w = build_rhs_mono_unstead_diff(operator_w, fw, capacity_w, bc[3], ic[3], vcat(u0z[1], u0z[2]), Δt, t)

    return b_u, b_v, b_w
end

function solve_DiffusionVecUnsteadyMono!(solver::SolverVec, phase::VectorPhase{1}, u0, Δt::Float64, Tend::Float64, bc::Tuple{BorderConditions}, ic::Tuple{AbstractBoundary}; method::Function = gmres, kwargs...)
    println("Résolution du problème:")
    println("- Vector problem")
    println("- Monophasic problem")
    println("- Unsteady problem")
    println("- Diffusion problem")
    
    # Initialisation
    t = 0.0
    u = vcat(u0[1], u0[2])
    n = 0

    lenu = length(u0[1])
    n_recon = Int(size(solver.A, 1) / 2)
    # Time loop
    while t < Tend
        n += 1
        t += Δt

        # Build the matrix and the right-hand side
        solver.A = build_mono_unstead_diff_vec_matrix(phase, bc, ic, Δt)
        solver.b = build_mono_unstead_diff_vec_rhs(phase, bc, ic, Δt, u0, t)

        if method == \
            # Remove zero rows and columns for direct solver
            solver.A, solver.b, rows_idx, cols_idx = remove_zero_rows_cols!(solver.A, solver.b)
            x_reduced = solver.A \ solver.b
            # Reconstruct the full solution vector
            solver.x = zeros(2*n_recon)
            solver.x[cols_idx] = x_reduced
        else
            # Use iterative solver directly
            kwargs_nt = (; kwargs...)
            log = get(kwargs_nt, :log, false)
            if log
                solver.x, solver.ch = method(solver.A, solver.b; kwargs...)
            else
                solver.x = method(solver.A, solver.b; kwargs...)
            end
        end
        # Update the solution
        push!(solver.states, solver.x)
        @show maximum(solver.x)
        u0 = (solver.x[1:lenu], solver.x[lenu+1:end])
    end

    return u
end

function solve_DiffusionVecUnsteadyMono!(solver::SolverVec, phase::VectorPhase{2}, u0x, u0y, Δt::Float64, Tend::Float64, bc::Tuple{BorderConditions, BorderConditions}, ic::Tuple{AbstractBoundary, AbstractBoundary}; method::Function = gmres, kwargs...)
    println("Résolution du problème:")
    println("- Vector problem")
    println("- Monophasic problem")
    println("- Unsteady problem")
    println("- Diffusion problem")
    
    # Initialisation
    t = 0.0
    u = vcat(u0x[1], u0x[2], u0y[1], u0y[2])
    n = 0

    lenu = length(u0x[1])
    lenv = length(u0y[1])

    n_recon = Int(size(solver.A, 1) / 4)
    # Time loop
    while t < Tend
        n += 1
        t += Δt

        # Build the matrix and the right-hand side
        A_u, A_v = build_mono_unstead_diff_vec_matrix(phase, bc, ic, Δt)
        b_u, b_v = build_mono_unstead_diff_vec_rhs(phase, bc, ic, Δt, u0x, u0y, t)

        solver.A = blockdiag(A_u, A_v)
        solver.b = vcat(b_u, b_v)

        if method == \
            # Remove zero rows and columns for direct solver
            solver.A, solver.b, rows_idx, cols_idx = remove_zero_rows_cols!(solver.A, solver.b)
            x_reduced = solver.A \ solver.b
            # Reconstruct the full solution vector
            solver.x = zeros(4*n_recon)
            solver.x[cols_idx] = x_reduced
        else
            # Use iterative solver directly
            kwargs_nt = (; kwargs...)
            log = get(kwargs_nt, :log, false)
            if log
                solver.x, solver.ch = method(solver.A, solver.b; kwargs...)
            else
                solver.x = method(solver.A, solver.b; kwargs...)
            end
        end

        # Update the solution
        push!(solver.states, solver.x)
        @show maximum(solver.x)
        u0x = (solver.x[1:lenu], solver.x[lenu+1:2*lenu])
        u0y = (solver.x[2*lenu+1:2*lenu+lenv], solver.x[2*lenu+lenv+1:end]) 
    end

    return u
end

function solve_DiffusionVecUnsteadyMono!(solver::SolverVec, phase::VectorPhase{3}, u0x, u0y, u0z, Δt::Float64, Tend::Float64, bc::Tuple{BorderConditions, BorderConditions, BorderConditions}, ic::Tuple{AbstractBoundary, AbstractBoundary, AbstractBoundary}; method::Function = gmres, kwargs...)
    println("Résolution du problème:")
    println("- Vector problem")
    println("- Monophasic problem")
    println("- Unsteady problem")
    println("- Diffusion problem")
    
    # Initialisation
    t = 0.0
    u = vcat(u0x[1], u0x[2], u0y[1], u0y[2], u0z[1], u0z[2])
    n = 0

    lenu = length(u0x[1])
    lenv = length(u0y[1])
    lenw = length(u0z[1])

    n_recon = Int(div(size(solver.A, 1),6))

    # Time loop
    while t < Tend
        n += 1
        t += Δt

        # Build the matrix and the right-hand side
        A_u, A_v, A_w = build_mono_unstead_diff_vec_matrix(phase, bc, ic, Δt)
        b_u, b_v, b_w = build_mono_unstead_diff_vec_rhs(phase, bc, ic, Δt, u0x, u0y, u0z, t)

        solver.A = blockdiag(A_u, A_v, A_w)
        solver.b = vcat(b_u, b_v, b_w)

        if method == \
            # Remove zero rows and columns for direct solver
            solver.A, solver.b, rows_idx, cols_idx = remove_zero_rows_cols!(solver.A, solver.b)
            x_reduced = solver.A \ solver.b
            # Reconstruct the full solution vector
            solver.x = zeros(6*n_recon)
            solver.x[cols_idx] = x_reduced
        else
            # Use iterative solver directly
            kwargs_nt = (; kwargs...)
            log = get(kwargs_nt, :log, false)
            if log
                solver.x, solver.ch = method(solver.A, solver.b; kwargs...)
            else
                solver.x = method(solver.A, solver.b; kwargs...)
            end
        end

        # Update the solution
        push!(solver.states, solver.x)
        @show maximum(solver.x)
        u0x = (solver.x[1:lenu], solver.x[lenu+1:2*lenu])
        u0y = (solver.x[2*lenu+1:2*lenu+lenv], solver.x[2*lenu+lenv+1:2*lenu+2*lenv])
        u0z = (solver.x[2*lenu+2*lenv+1:2*lenu+2*lenv+lenw], solver.x[2*lenu+2*lenv + lenw+1:end])
    end

    return u
end






























# Vector Convection - Unsteady - Monophasic
function ConvectionVecUnsteadyMono(phase::VectorPhase{2}, bc::Tuple{BorderConditions, BorderConditions}, ic::Tuple{AbstractBoundary, AbstractBoundary}, Δt::Float64, Tend::Float64, u0::Vector{Float64})
    println("Création du solveur:")
    println("- Vector problem")
    println("- Monophasic problem")
    println("- Unsteady problem")
    println("- Convection problem")
    
    s = SolverVec(Unsteady, Monophasic, Advection, nothing, nothing, nothing, ConvergenceHistory(), [])

    A_u, A_v = build_mono_unstead_conv_vec_matrix(phase, bc, ic, Δt)
    b_u, b_v = build_mono_unstead_conv_vec_rhs(phase, bc, ic, Δt, u0, 0.0)

    s.A = blockdiag(A_u, A_v)
    s.b = vcat(b_u, b_v)
    # BC_border_mono!(s, phase, bc)

    return s
end

function build_mono_unstead_conv_vec_matrix(phase, bc, ic, Δt)
    operator_u, operator_v = phase.operator
    capacity_u, capacity_v = phase.capacity

    # For u-component
    #A_u = build_mono_unstead_adv_matrix(operator_u, capacity_u, bc[1], ic[1], Δt)
    A_u = build_mono_unstead_adv_diff_matrix(operator_u, capacity_u, 1.0, bc[1], ic[1], Δt)

    # For v-component
    #A_v = build_mono_unstead_adv_matrix(operator_v, capacity_v, bc[2], ic[2], Δt)
    A_v = build_mono_unstead_adv_diff_matrix(operator_v, capacity_v, 1.0, bc[2], ic[2], Δt)

    return A_u, A_v
end

function build_mono_unstead_conv_vec_rhs(phase, bc, ic, Δt, u0, t)
    operator_u, operator_v = phase.operator
    capacity_u, capacity_v = phase.capacity
    fu, fv = phase.source

    # For u-component
    #b_u = build_rhs_mono_unstead_adv(operator_u, fu, capacity_u, bc[1], ic[1], u0[1:end÷2], Δt, t)
    b_u = build_rhs_mono_unstead_adv_diff(operator_u, fu, capacity_u, bc[1], ic[1], u0[1:end÷2], Δt, t)

    # For v-component
    #b_v = build_rhs_mono_unstead_adv(operator_v, fv, capacity_v, bc[2], ic[2], u0[end÷2+1:end], Δt, t)
    b_v = build_rhs_mono_unstead_adv_diff(operator_v, fv, capacity_v, bc[2], ic[2], u0[end÷2+1:end], Δt, t)

    return b_u, b_v
end

function solve_ConvectionVecUnsteadyMono!(solver::SolverVec, phase::VectorPhase, u0::Vector{Float64}, Δt::Float64, Tend::Float64, bc::Tuple{BorderConditions, BorderConditions}, ic::Tuple{AbstractBoundary, AbstractBoundary}; method=IterativeSolvers.gmres, kwargs...)
    println("Résolution du problème:")
    println("- Vector problem")
    println("- Monophasic problem")
    println("- Unsteady problem")
    println("- Convection problem")
    
    # Initialisation
    t = 0.0
    u = copy(u0)
    n = 0
    n_recon = Int(size(solver.A, 1) / 4)

    # Time loop
    while t < Tend
        n += 1
        t += Δt

        # Build the matrix and the right-hand side
        A_u, A_v = build_mono_unstead_conv_vec_matrix(phase, bc, ic, Δt)
        b_u, b_v = build_mono_unstead_conv_vec_rhs(phase, bc, ic, Δt, u, t)

        solver.A = blockdiag(A_u, A_v)
        solver.b = vcat(b_u, b_v)

        # Solve the linear system
        if method == \
            # Remove zero rows and columns for direct solver
            solver.A, solver.b, rows_idx, cols_idx = remove_zero_rows_cols!(solver.A, solver.b)
            x_reduced = solver.A \ solver.b
            # Reconstruct the full solution vector
            solver.x = zeros(4*n_recon)
            solver.x[cols_idx] = x_reduced
        else
            # Use iterative solver directly
            kwargs_nt = (; kwargs...)
            log = get(kwargs_nt, :log, false)
            if log
                solver.x, solver.ch = method(solver.A, solver.b; kwargs...)
            else
                solver.x = method(solver.A, solver.b; kwargs...)
            end
        end

        # Update the solution
        push!(solver.states, solver.x)
        @show maximum(solver.x)
        u = copy(solver.x)
    end

    return u
end