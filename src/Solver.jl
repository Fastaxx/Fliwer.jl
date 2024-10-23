@enum TimeType begin
    Steady  # ∂ₜT = 0
    Unsteady # ∂ₜT ≠ 0
end

@enum PhaseType begin
    Monophasic  # Single phase
    Diphasic    # Two phases
end

@enum EquationType begin
    Diffusion           # ∂ₜT = ∇·(∇T) + S
    Advection           # ∂ₜT = -∇·(uT) + S
    DiffusionAdvection  # ∂ₜT = ∇·(D∇T) - ∇·(uT) + S
end

struct Solver{TT<:TimeType, PT<:PhaseType, ET<:EquationType}
    time_type::TT
    phase_type::PT
    equation_type::ET
    diffusion_coefficient::Float64  # Rajouter la possibilité d'une fonction
    advection_coefficient::Float64  # Rajouter la possibilité d'une fonction
end

function DiffusionSolver(steady::Bool, monophasic::Bool, D::Float64)
    time_type = steady ? Steady : Unsteady
    phase_type = monophasic ? Monophasic : Diphasic
    equation_type = Diffusion
    return Solver(time_type, phase_type, equation_type, D, 0.0) 
end
        
function AdvectionSolver(steady::Bool, monophasic::Bool, u::Float64)
    time_type = steady ? Steady : Unsteady
    phase_type = monophasic ? Monophasic : Diphasic
    equation_type = Advection
    return Solver(time_type, phase_type, equation_type, 0.0, u)
end

function AdvectionDiffusionSolver(steady::Bool, monophasic::Bool, D::Float64, u::Float64)
    time_type = steady ? Steady : Unsteady
    phase_type = monophasic ? Monophasic : Diphasic
    equation_type = DiffusionAdvection
    return Solver(time_type, phase_type, equation_type, D, u)
end

function initialize(s::Solver)
    println("Initialisation du solveur:")

    # EquationType
    if s.equation_type == Diffusion
        println("  - Mode: Diffusion uniquement")
        # Initialiser les paramètres spécifiques à la diffusion
    elseif s.equation_type == Advection
        println("  - Mode: Advection uniquement")
        # Initialiser les paramètres spécifiques à l'advection
    elseif s.equation_type == DiffusionAdvection
        println("  - Mode: Advection-Diffusion")
        # Initialiser les paramètres spécifiques à l'advection-diffusion
    end

    # TimeType
    if s.time_type == Steady
        println("  - Simulation : Steady")
        # Initialiser les paramètres spécifiques à une simulation stationnaire
    elseif s.time_type == Unsteady
        println("  - Simulation : Transient")
        # Initialiser les paramètres spécifiques à une simulation transitoire
    end

    # PhaseType
    if s.phase_type == Monophasic
        println("  - Phase : Monophasic")
        # Initialiser les paramètres spécifiques à une phase unique
    elseif s.phase_type == Diphasic
        println("  - Phase : Diphasic")
        # Initialiser les paramètres spécifiques à deux phases
    end

end
