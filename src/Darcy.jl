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
end
