# Monophasic Steady
function solve!(s::Solver, phase::Phase; method::Function = gmres, kwargs...)
    if s.A === nothing
        error("Solver is not initialized. Call a solver constructor first.")
    end

    kwargs_nt = (; kwargs...)
    log = get(kwargs_nt, :log, false)

    if log
        s.x, s.ch = method(s.A, s.b; kwargs...)
    else
        s.x = method(s.A, s.b; kwargs...)
    end
end

# Diphasic Steady
function solve!(s::Solver, phase1::Phase, phase2::Phase; method::Function = gmres, kwargs...)
    if s.A === nothing
        error("Solver is not initialized. Call a solver constructor first.")
    end

    kwargs_nt = (; kwargs...)
    log = get(kwargs_nt, :log, false)

    if log
        s.x, s.ch = method(s.A, s.b; kwargs...)
    else
        s.x = method(s.A, s.b; kwargs...)
    end
end

# Monophasic Unsteady
function solve!(s::Solver, phase::Phase, Tᵢ, Δt::Float64, Tₑ, bc_b::BorderConditions, bc::AbstractBoundary; method::Function = gmres, kwargs...)
    if s.A === nothing
        error("Solver is not initialized. Call a solver constructor first.")
    end

    kwargs_nt = (; kwargs...)
    log = get(kwargs_nt, :log, false)

    BC_border_mono!(s.A, s.b, bc_b, phase.capacity.mesh)
    if log
        s.x, s.ch = method(s.A, s.b; kwargs...)
    else
        s.x = method(s.A, s.b; kwargs...)
    end
    t = 0.0
    while t < Tₑ
        t += Δt
        println("Time: ", t)
        s.b = build_rhs(phase.operator, phase.source, phase.capacity, bc_b, bc, Tᵢ, Δt, t)
        BC_border_mono!(s.A, s.b, bc_b, phase.capacity.mesh)
        if log
            s.x, s.ch = method(s.A, s.b; kwargs...)
        else
            s.x = method(s.A, s.b; kwargs...)
        end
        push!(s.states, s.x)
        @show maximum(s.x)
        Tᵢ = s.x
    end
end

# Diphasic Unsteady
function solve!(s::Solver, phase1::Phase, phase2::Phase, Tᵢ, Δt::Float64, Tₑ, bc_b::BorderConditions, ic::InterfaceConditions; method::Function = gmres, kwargs...)
    if s.A === nothing
        error("Solver is not initialized. Call a solver constructor first.")
    end

    kwargs_nt = (; kwargs...)
    log = get(kwargs_nt, :log, false)

    BC_border_diph!(s.A, s.b, bc_b, phase2.capacity.mesh)
    if log
        s.x, s.ch = method(s.A, s.b; kwargs...)
    else
        s.x = method(s.A, s.b; kwargs...)
    end
    t = 0.0
    while t < Tₑ
        t += Δt
        println("Time: ", t)
        s.b = build_rhs(phase1.operator, phase2.operator, phase1.source, phase2.source, phase1.capacity, phase2.capacity, bc_b, ic, Tᵢ, Δt, t)
        BC_border_diph!(s.A, s.b, bc_b, phase2.capacity.mesh)
        if log
            s.x, s.ch = method(s.A, s.b; kwargs...)
        else
            s.x = method(s.A, s.b; kwargs...)
        end
        push!(s.states, s.x)
        Tᵢ = s.x
    end
end
