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

mutable struct Solver{TT<:TimeType, PT<:PhaseType, ET<:EquationType}
    time_type::TT
    phase_type::PT
    equation_type::ET
    diffusion_coefficient::Float64  # Rajouter la possibilité d'une fonction
    A::Union{SparseMatrixCSC{Float64, Int}, Nothing}
    b::Union{Vector{Float64}, Nothing}
end

function DiffusionSolver(steady::Bool, monophasic::Bool, D::Float64)
    time_type = steady ? Steady : Unsteady
    phase_type = monophasic ? Monophasic : Diphasic
    equation_type = Diffusion
    return Solver(time_type, phase_type, equation_type, D, nothing, nothing)
end

function initialize!(s::Solver, operator::AbstractOperators, bc::BorderCondition)
    println("Initialisation du solveur:")

    # EquationType
    if s.equation_type == Diffusion
        println("  - Mode: Diffusion uniquement")
        s.A = build_diffusion_matrix(operator, s.diffusion_coefficient, bc)
    elseif s.equation_type == Advection
        println("  - Mode: Advection uniquement")
        s.A = build_advection_matrix(operator, s.advection_coefficient, bc)
    elseif s.equation_type == DiffusionAdvection
        println("  - Mode: Advection-Diffusion")
        s.A = build_diffusion_advection_matrix(operator, s.diffusion_coefficient, s.advection_coefficient, bc)
    end

    # TimeType
    if s.time_type == Steady
        println("  - Simulation : Steady")
        # Pour les simulations steady, aucune initialisation supplémentaire
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

# Fonction pour construire la matrice de diffusion
function build_diffusion_matrix(operator::DiffusionOps, D::Float64, bc::Dirichlet)
    N = prod(operator.size)

    Iᵦ = spzeros(N, N)
    Iₐ = I(N)
    Iᵧ = spzeros(N, N)
    
    A = vcat(hcat(operator.G' * operator.Wꜝ * operator.G, operator.G' * operator.Wꜝ * operator.H), hcat(Iᵦ * operator.H' * operator.Wꜝ * operator.G, Iᵦ * operator.H' * operator.Wꜝ * operator.H + Iₐ * Iᵧ))

    #BC_border!(A, bc)
    #BC_interface!(A, bc)
    return A
end

function build_rhs(operator::DiffusionOps, bc, f)
    N = prod(operator.size)
    fₒ = zeros(N)
    gᵧ = ones(N)
    b = zeros(2N)

    Iᵧ = spzeros(N, N)

    # Compute the source term
    for i in 1:N
        x, y, z = 0., 0., 0.
        fₒ[i] = f(x, y, z)
    end
    
    b = vcat(operator.V*fₒ, Iᵧ * gᵧ)

    #BC_border_b!(b, bc)
    #BC_interface_b!(b, bc)

    return b
end


function solve!(s::Solver, operator::DiffusionOps, bc, f)
    if s.A === nothing
        error("Solver is not initialized. Call initialize!(solver, operator, bc) first.")
    end

    b = build_rhs(operator, bc, f)

    s.b = b

    T = cg(s.A, b)

    return T
end
