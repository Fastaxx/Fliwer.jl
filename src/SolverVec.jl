# Vector Diffusion - Unsteady - Monophasic
function DiffusionVecUnsteadyMono(phase::VectorPhase, bc::Tuple{BorderConditions, BorderConditions}, ic::Tuple{AbstractBoundary, AbstractBoundary}, Δt::Float64, Tend::Float64, u0::Vector{Float64})
    println("Création du solveur:")
    println("- Vector problem")
    println("- Monophasic problem")
    println("- Unsteady problem")
    println("- Diffusion problem")
    
    s = Solver(Unsteady, Monophasic, Diffusion, nothing, nothing, nothing, ConvergenceHistory(), [])

    A_u, A_v = build_mono_unstead_diff_vec_matrix(phase, bc, ic, Δt)
    b_u, b_v = build_mono_unstead_diff_vec_rhs(phase, bc, ic, Δt, u0, 0.0)

    s.A = blockdiag(A_u, A_v)
    s.b = vcat(b_u, b_v)
    # BC_border_mono!(s, phase, bc)

    return s
end

function build_mono_unstead_diff_vec_matrix(phase, bc, ic, Δt)
    operator_u, operator_v = phase.operator
    capacity_u, capacity_v = phase.capacity
    μ = phase.Diffusion_coeff

    # For u-component
    A_u = build_mono_unstead_diff_matrix(operator_u, capacity_u, μ, bc[1], ic[1], Δt)

    # For v-component
    A_v = build_mono_unstead_diff_matrix(operator_v, capacity_v, μ, bc[2], ic[2], Δt)

    return A_u, A_v
end

function build_mono_unstead_diff_vec_rhs(phase, bc, ic, Δt, u0, t)
    operator_u, operator_v = phase.operator
    capacity_u, capacity_v = phase.capacity
    fu, fv = phase.source
    μ = phase.Diffusion_coeff

    # For u-component
    b_u = build_rhs_mono_unstead_diff(operator_u, fu, capacity_u, bc[1], ic[1], u0[1:end÷2], Δt, t)

    # For v-component
    b_v = build_rhs_mono_unstead_diff(operator_v, fv, capacity_v, bc[2], ic[2], u0[end÷2+1:end], Δt, t)

    return b_u, b_v
end

function solve_DiffusionVecUnsteadyMono!(solver::Solver, phase::VectorPhase, u0::Vector{Float64}, Δt::Float64, Tend::Float64, bc::Tuple{BorderConditions, BorderConditions}, ic::Tuple{AbstractBoundary, AbstractBoundary}; method=IterativeSolvers.bicgstabl)
    println("Résolution du problème:")
    println("- Vector problem")
    println("- Monophasic problem")
    println("- Unsteady problem")
    println("- Diffusion problem")
    
    # Initialisation
    t = 0.0
    u = copy(u0)
    n = 0

    # Time loop
    while t < Tend
        n += 1
        t += Δt

        # Build the matrix and the right-hand side
        A_u, A_v = build_mono_unstead_diff_vec_matrix(phase, bc, ic, Δt)
        b_u, b_v = build_mono_unstead_diff_vec_rhs(phase, bc, ic, Δt, u, t)

        solver.A = blockdiag(A_u, A_v)
        solver.b = vcat(b_u, b_v)

        # Solve the linear system
        solver.x = method(solver.A, solver.b)

        # Update the solution
        push!(solver.states, solver.x)
        @show maximum(solver.x)
        u = copy(solver.x)
    end

    return u
end