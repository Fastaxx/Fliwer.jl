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
    A::Union{SparseMatrixCSC{Float64, Int}, Nothing}
    b::Union{Vector{Float64}, Nothing}
    x::Union{Vector{Float64}, Nothing}
    states::Vector{Any}
end


# Constructeur de solveur : Diffusion - Steady - Monophasic
function DiffusionSteadyMono(phase::Phase, bc_b::BorderConditions, bc_i::AbstractBoundary)
    println("Création du solveur:")
    println("- Monophasic problem")
    println("- Steady problem")
    println("- Diffusion problem")
    
    s = Solver(Steady, Monophasic, Diffusion, nothing, nothing, nothing, [])
    
    s.A = build_mono_stead_diff_matrix(phase.operator, phase.capacity, phase.Diffusion_coeff, bc_b, bc_i)
    s.b = build_rhs(phase.operator, phase.source, phase.capacity, bc_b, bc_i)

    return s
end

function build_mono_stead_diff_matrix(operator::DiffusionOps, capacity::Capacity, D::Float64, bc_b::BorderConditions, bc::AbstractBoundary)
    n = prod(operator.size)
    Iₐ, Iᵦ = build_I_bc(operator, bc)
    Iᵧ = build_I_g(operator)
    Id = build_I_D(operator, D)

    A = vcat(hcat(Id * operator.G' * operator.Wꜝ * operator.G, Id * operator.G' * operator.Wꜝ * operator.H), hcat(Iᵦ * operator.H' * operator.Wꜝ * operator.G, Iᵦ * operator.H' * operator.Wꜝ * operator.H + Iₐ * Iᵧ))

    #BC_border!(A, bc_b, capacity.mesh.tag)
    return A
end

function build_rhs(operator::DiffusionOps, f, capacite::Capacity, bc_b::BorderConditions, bc::AbstractBoundary)
    N = prod(operator.size)
    b = zeros(2N)

    Iᵧ = build_I_g(operator)
    fₒ = build_source(operator, f, capacite)
    gᵧ = build_g_g(operator, bc, capacite)

    # Build the right-hand side
    b = vcat(operator.V*fₒ, Iᵧ * gᵧ)

    #BC_border_b!(b, bc_b)
    #BC_interface_b!(b, bc)

    return b
end

function solve!(s::Solver, phase::Phase)
    if s.A === nothing
        error("Solver is not initialized. Call a solver constructor first.")
    end

    s.x = gmres(s.A, s.b, abstol=1e-15)
end


# Constructeur de solveur : Diffusion - Steady - Diphasic
function DiffusionSteadyDiph(phase1::Phase, phase2::Phase, bc_b::BorderConditions, ic::InterfaceConditions)
    println("Création du solveur:")
    println("- Diphasic problem")
    println("- Steady problem")
    println("- Diffusion problem")
    
    s = Solver(Steady, Diphasic, Diffusion, nothing, nothing, nothing, [])
    
    s.A = build_diph_stead_diff_matrix(phase1.operator, phase2.operator, phase1.Diffusion_coeff, phase2.Diffusion_coeff, bc_b, ic)
    s.b = build_rhs(phase1.operator, phase2.operator, phase1.source, phase2.source, phase1.capacity, phase2.capacity, bc_b, ic)

    return s
end

function build_diph_stead_diff_matrix(operator1::DiffusionOps, operator2::DiffusionOps, D1::Float64, D2::Float64, bc_b::BorderConditions, ic::InterfaceConditions)
    n = prod(operator1.size)

    jump, flux = ic.scalar, ic.flux
    Iₐ1, Iₐ2 = jump.α₁*I(n), jump.α₂*I(n)
    Iᵦ1, Iᵦ2 = flux.β₁*I(n), flux.β₂*I(n)
    Id1, Id2 = build_I_D(operator1, D1), build_I_D(operator2, D2)


    block1 = Id1 * operator1.G' * operator1.Wꜝ * operator1.G
    block2 = Id1 * operator1.G' * operator1.Wꜝ * operator1.H
    block3 = Id2 * operator2.G' * operator2.Wꜝ * operator2.G
    block4 = Id2 * operator2.G' * operator2.Wꜝ * operator2.H
    block5 = operator1.H' * operator1.Wꜝ * operator1.G
    block6 = operator1.H' * operator1.Wꜝ * operator1.H 
    block7 = operator2.H' * operator2.Wꜝ * operator2.G
    block8 = operator2.H' * operator2.Wꜝ * operator2.H

    A = vcat(hcat(block1, block2, zeros(n, n), zeros(n, n)),
             hcat(zeros(n, n), Iₐ1, zeros(n, n), -Iₐ2),
             hcat(zeros(n, n), zeros(n, n), block3, block4),
             hcat(Iᵦ1*block5, Iᵦ1*block6, Iᵦ2*block7, Iᵦ2*block8))

    #BC_border!(A, bc_b1, bc_b2)
    return A
end

function build_rhs(operator1::DiffusionOps, operator2::DiffusionOps, f1, f2, capacite1::Capacity, capacite2::Capacity, bc_b::BorderConditions, ic::InterfaceConditions)
    N = prod(operator1.size)
    b = zeros(4N)

    jump, flux = ic.scalar, ic.flux
    Iᵧ1, Iᵧ2 = build_I_g(operator1), build_I_g(operator2)
    gᵧ, hᵧ = build_g_g(operator1, jump, capacite1), build_g_g(operator2, flux, capacite2)

    fₒ1 = build_source(operator1, f1, capacite1)
    fₒ2 = build_source(operator2, f2, capacite2)

    # Build the right-hand side
    b = vcat(operator1.V*fₒ1, gᵧ, operator2.V*fₒ2, Iᵧ2*hᵧ)

    #BC_border_b!(b, bc_b)
    #BC_interface_b!(b, ic)

    return b
end

function solve!(s::Solver, phase1::Phase, phase2::Phase)
    if s.A === nothing
        error("Solver is not initialized. Call a solver constructor first.")
    end

    s.x = bicgstabl(s.A, s.b, abstol=1e-15)
end


# Constructeur de solveur : Diffusion - Unsteady - Monophasic
function DiffusionUnsteadyMono(phase::Phase, bc_b::BorderConditions, bc_i::AbstractBoundary, Δt::Float64, Tₑ::Float64, Tᵢ::Vector{Float64})
    println("Création du solveur:")
    println("- Monophasic problem")
    println("- Unsteady problem")
    println("- Diffusion problem")
    
    s = Solver(Unsteady, Monophasic, Diffusion, nothing, nothing, nothing, [])
    
    s.A = build_mono_unstead_diff_matrix(phase.operator, phase.Diffusion_coeff, bc_b, bc_i, Δt)
    s.b = build_rhs(phase.operator, phase.source, phase.capacity, bc_b, bc_i, Tᵢ, Δt, 0.0)

    return s
end

function build_mono_unstead_diff_matrix(operator::DiffusionOps, D::Float64, bc_b::BorderConditions, bc::AbstractBoundary, Δt::Float64)
    n = prod(operator.size)
    Iₐ, Iᵦ = build_I_bc(operator, bc)
    Iᵧ = build_I_g(operator)
    Id = build_I_D(operator, D)

    block1 = operator.V + Δt/2 * Id * operator.G' * operator.Wꜝ * operator.G
    block2 = Δt/2 * Id * operator.G' * operator.Wꜝ * operator.H
    block3 = Iᵦ * operator.H' * operator.Wꜝ * operator.G
    block4 = Iᵦ * operator.H' * operator.Wꜝ * operator.H + Iₐ * Iᵧ

    A = vcat(hcat(block1, block2), hcat(block3, block4))

    #BC_border!(A, bc_b)
    return A
end

function build_rhs(operator::DiffusionOps, f, capacite::Capacity, bc_b::BorderConditions, bc::AbstractBoundary, Tᵢ, Δt::Float64, t::Float64)
    N = prod(operator.size)
    b = zeros(2N)

    Iᵧ = build_I_g(operator)
    fₒn, fₒn1 = build_source(operator, f, t, capacite), build_source(operator, f, t+Δt, capacite)
    gᵧ = build_g_g(operator, bc, capacite)

    Tₒ, Tᵧ = Tᵢ[1:N], Tᵢ[N+1:end]

    # Build the right-hand side
    b = vcat((operator.V - Δt/2 * operator.G' * operator.Wꜝ * operator.G)*Tₒ - Δt/2 * operator.G' * operator.Wꜝ * operator.H * Tᵧ + Δt/2 * operator.V * (fₒn + fₒn1), Iᵧ * gᵧ)

    #BC_border_b!(b, bc_b)
    #BC_interface_b!(b, bc)

    return b
end

function solve!(s::Solver, phase::Phase, Tᵢ, Δt::Float64, Tₑ, bc_b::BorderConditions, bc::AbstractBoundary)
    if s.A === nothing
        error("Solver is not initialized. Call a solver constructor first.")
    end

    s.x = cg(s.A, s.b)
    t=0.0
    while t < Tₑ
        t+=Δt
        println("Time: ", t)
        s.b = build_rhs(phase.operator, phase.source, phase.capacity, bc_b, bc, Tᵢ, Δt, t)
        
        s.x = cg(s.A, s.b)
        push!(s.states, s.x)
        @show maximum(s.x)

        Tᵢ = s.x

    end
end


# Constructeur de solveur : Diffusion - Unsteady - Diphasic
function DiffusionUnsteadyDiph(phase1::Phase, phase2::Phase, bc_b::BorderConditions, ic::InterfaceConditions, Δt::Float64, Tₑ::Float64, Tᵢ::Vector{Float64})
    println("Création du solveur:")
    println("- Diphasic problem")
    println("- Unsteady problem")
    println("- Diffusion problem")
    
    s = Solver(Unsteady, Diphasic, Diffusion, nothing, nothing, nothing, [])
    
    s.A = build_diph_unstead_diff_matrix(phase1.operator, phase2.operator, phase1.Diffusion_coeff, phase2.Diffusion_coeff, bc_b, ic, Δt)
    s.b = build_rhs(phase1.operator, phase2.operator, phase1.source, phase2.source, phase1.capacity, phase2.capacity, bc_b, ic, Tᵢ, Δt, 0.0)

    return s
end

function build_diph_unstead_diff_matrix(operator1::DiffusionOps, operator2::DiffusionOps, D1::Float64, D2::Float64, bc_b::BorderConditions, ic::InterfaceConditions, Δt::Float64)
    n = prod(operator1.size)

    jump, flux = ic.scalar, ic.flux
    Iₐ1, Iₐ2 = jump.α₁*I(n), jump.α₂*I(n)
    Iᵦ1, Iᵦ2 = flux.β₁*I(n), flux.β₂*I(n)
    Id1, Id2 = build_I_D(operator1, D1), build_I_D(operator2, D2)

    block1 = operator1.V + Δt/2 * Id1 * operator1.G' * operator1.Wꜝ * operator1.G
    block2 = Δt/2 * Id1 * operator1.G' * operator1.Wꜝ * operator1.H
    block3 = operator2.V + Δt/2 * Id2 * operator2.G' * operator2.Wꜝ * operator2.G
    block4 = Δt/2 * Id2 * operator2.G' * operator2.Wꜝ * operator2.H
    block5 = Iᵦ1 * operator1.H' * operator1.Wꜝ * operator1.G
    block6 = Iᵦ1 * operator1.H' * operator1.Wꜝ * operator1.H
    block7 = Iᵦ2 * operator2.H' * operator2.Wꜝ * operator2.G
    block8 = Iᵦ2 * operator2.H' * operator2.Wꜝ * operator2.H

    A = vcat(hcat(block1, block2, zeros(n, n), zeros(n, n)),
             hcat(zeros(n, n), Iₐ1, zeros(n, n), -Iₐ2),
             hcat(zeros(n, n), zeros(n, n), block3, block4),
             hcat(Iᵦ1*block5, Iᵦ1*block6, Iᵦ2*block7, Iᵦ2*block8))

    #BC_border!(A, bc_b1, bc_b2)
    return A
end

function build_rhs(operator1::DiffusionOps, operator2::DiffusionOps, f1, f2, capacite1::Capacity, capacite2::Capacity, bc_b::BorderConditions, ic::InterfaceConditions, Tᵢ, Δt::Float64, t::Float64)
    N = prod(operator1.size)
    b = zeros(4N)

    jump, flux = ic.scalar, ic.flux
    Iᵧ1, Iᵧ2 = build_I_g(operator1), build_I_g(operator2)
    gᵧ, hᵧ = build_g_g(operator1, jump,capacite1), build_g_g(operator2, flux, capacite2)

    fₒn1, fₒn2 = build_source(operator1, f1, t, capacite1), build_source(operator2, f2, t, capacite2)
    fₒn1p1, fₒn2p1 = build_source(operator1, f1, t+Δt, capacite1), build_source(operator2, f2, t+Δt, capacite2)

    Tₒ1, Tᵧ1 = Tᵢ[1:N], Tᵢ[N+1:2N]
    Tₒ2, Tᵧ2 = Tᵢ[2N+1:3N], Tᵢ[3N+1:end]

    # Build the right-hand side
    b = vcat((operator1.V - Δt/2 * operator1.G' * operator1.Wꜝ * operator1.G)*Tₒ1 - Δt/2 * operator1.G' * operator1.Wꜝ * operator1.H * Tᵧ1 + Δt/2 * operator1.V * (fₒn1 + fₒn1p1), gᵧ, (operator2.V - Δt/2 * operator2.G' * operator2.Wꜝ * operator2.G)*Tₒ2 - Δt/2 * operator2.G' * operator2.Wꜝ * operator2.H * Tᵧ2 + Δt/2 * operator2.V * (fₒn2 + fₒn2p1), Iᵧ2*hᵧ)

    #BC_border_b!(b, bc_b)
    #BC_interface_b!(b, ic)

    return b
end

function solve!(s::Solver, phase1::Phase, phase2::Phase, Tᵢ, Δt::Float64, Tₑ, bc_b::BorderConditions, ic::InterfaceConditions)
    if s.A === nothing
        error("Solver is not initialized. Call a solver constructor first.")
    end

    s.x = cg(s.A, s.b)
    t=0.0
    while t < Tₑ
        t+=Δt
        println("Time: ", t)
        s.b = build_rhs(phase1.operator, phase2.operator, phase1.source, phase2.source, phase1.capacity, phase2.capacity, bc_b, ic, Tᵢ, Δt, t)
        
        s.x = gmres(s.A, s.b)
        push!(s.states, s.x)
        

        Tᵢ = s.x

    end
end







    
function build_I_D(operator::DiffusionOps, D::Union{Float64,Function})
    n = prod(operator.size)
    Id = spdiagm(0 => ones(n))

    if D isa Function
        Id = D(Id)
    else
        Id = D * Id
    end
    return Id
end

function build_I_bc(operator,bc::AbstractBoundary)
    n = prod(operator.size)
    Iᵦ = spzeros(n, n)
    Iₐ = spzeros(n, n)

    if bc isa Dirichlet
        Iₐ = I(n)
    elseif bc isa Neumann
        Iᵦ = I(n)
    elseif bc isa Robin
        if bc.α isa Function
            Iₐ = bc.α(I(n))
        else
            Iₐ = bc.α * I(n)
            Iᵦ = bc.β * I(n)
        end 
    end
    return Iₐ, Iᵦ
end

function build_I_g(operator::AbstractOperators)
    n = prod(operator.size)
    vec_1 = [1 for i in 1:size(operator.H', 2)]
    row_sums = operator.H'*vec_1
    abs_row_sums = abs.(row_sums)
    #abs_row_sums[abs_row_sums .< 10^-8] .= 0.0  # Mettre à zéro les petites valeurs
    Iᵧ = spdiagm(0 => abs_row_sums)
    return Iᵧ
end

function build_g_g(operator::DiffusionOps, bc::Union{AbstractBoundary, AbstractInterfaceBC}, capacite::Capacity)
    n = prod(operator.size)
    gᵧ = ones(n)

    if bc.value isa Function
        for i in 1:n
            x, y, z = get_coordinates(i, capacite.C_γ)
            gᵧ[i] = bc.value(x, y, z)
        end
    else
        gᵧ = bc.value * gᵧ
    end
    return gᵧ
end

function build_source(operator::DiffusionOps, f, capacite::Capacity)
    N = prod(operator.size)
    fₒ = zeros(N)

    # Compute the source term
    for i in 1:N
        x, y, z = get_coordinates(i, capacite.C_ω)
        fₒ[i] = f(x, y, z)
    end

    return fₒ
end

function build_source(operator::DiffusionOps, f, t, capacite::Capacity)
    N = prod(operator.size)
    fₒ = zeros(N)

    # Compute the source term
    for i in 1:N
        x, y, z = get_coordinates(i, capacite.C_ω)
        fₒ[i] = f(x, y, z, t)
    end

    return fₒ
end

function get_coordinates(i, C_ω)
    if length(C_ω[1]) == 1
        x = C_ω[i][1]
        return x, 0., 0.
    elseif length(C_ω[1]) == 2
        x, y = C_ω[i][1], C_ω[i][2]
        return x, y, 0.
    else
        x, y, z = C_ω[i][1], C_ω[i][2], C_ω[i][3]
        return x, y, z
    end
end
