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

"""
    mutable struct Solver{TT<:TimeType, PT<:PhaseType, ET<:EquationType}

The `Solver` struct represents a solver for a specific type of problem.

# Fields
- `time_type::TT`: The type of time used in the solver : `Steady` or `Unsteady`.
- `phase_type::PT`: The type of phase used in the solver : `Monophasic` or `Diphasic`.
- `equation_type::ET`: The type of equation used in the solver : `Diffusion`, `Advection` or `DiffusionAdvection`.
- `A::Union{SparseMatrixCSC{Float64, Int}, Nothing}`: The coefficient matrix A of the equation system, if applicable.
- `b::Union{Vector{Float64}, Nothing}`: The right-hand side vector b of the equation system, if applicable.
- `x::Union{Vector{Float64}, Nothing}`: The solution vector x of the equation system, if applicable.
- `states::Vector{Any}`: The states of the system at different times, if applicable.

"""
mutable struct Solver{TT<:TimeType, PT<:PhaseType, ET<:EquationType}
    time_type::TT
    phase_type::PT
    equation_type::ET
    A::Union{SparseMatrixCSC{Float64, Int}, Nothing}
    b::Union{Vector{Float64}, Nothing}
    x::Union{Vector{Float64}, Nothing}
    ch::IterativeSolvers.ConvergenceHistory
    states::Vector{Any}
end


# Diffusion - Steady - Monophasic
"""
    DiffusionSteadyMono(phase::Phase, bc_b::BorderConditions, bc_i::AbstractBoundary)

Create a solver for a steady-state monophasic diffusion problem.

Arguments:
- `phase` : Phase object representing the phase of the problem.
- `bc_b` : BorderConditions object representing the boundary conditions of the problem at the boundary.
- `bc_i` : AbstractBoundary object representing the internal boundary conditions of the problem.

Returns:
- `s` : Solver object representing the created solver.
"""
function DiffusionSteadyMono(phase::Phase, bc_b::BorderConditions, bc_i::AbstractBoundary)
    println("Création du solveur:")
    println("- Monophasic problem")
    println("- Steady problem")
    println("- Diffusion problem")
    
    s = Solver(Steady, Monophasic, Diffusion, nothing, nothing, nothing, ConvergenceHistory(), [])
    
    s.A = build_mono_stead_diff_matrix(phase.operator, phase.capacity, phase.Diffusion_coeff, bc_b, bc_i)
    s.b = build_rhs_mono_stead_diff(phase.operator, phase.source, phase.capacity, bc_b, bc_i)

    BC_border_mono!(s.A, s.b, bc_b, phase.capacity.mesh)

    return s
end

function build_mono_stead_diff_matrix(operator::DiffusionOps, capacity::Capacity, D::Float64, bc_b::BorderConditions, bc::AbstractBoundary)
    n = prod(operator.size)
    Iₐ, Iᵦ = build_I_bc(operator, bc)
    Iᵧ =  build_I_g(operator) #capacity.Γ
    Id = build_I_D(operator, D)

    A = vcat(hcat(Id * operator.G' * operator.Wꜝ * operator.G, Id * operator.G' * operator.Wꜝ * operator.H), hcat(Iᵦ * operator.H' * operator.Wꜝ * operator.G, Iᵦ * operator.H' * operator.Wꜝ * operator.H + Iₐ * Iᵧ))

    return A
end

function build_rhs_mono_stead_diff(operator::DiffusionOps, f, capacite::Capacity, bc_b::BorderConditions, bc::AbstractBoundary)
    N = prod(operator.size)
    b = zeros(2N)

    Iᵧ = build_I_g(operator) #capacite.Γ 
    fₒ = build_source(operator, f, capacite)
    gᵧ = build_g_g(operator, bc, capacite)

    # Build the right-hand side
    b = vcat(operator.V*fₒ, Iᵧ * gᵧ)

    return b
end

function solve_DiffusionSteadyMono!(s::Solver, phase::Phase; method::Function = gmres, kwargs...)
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



# Diffusion - Steady - Diphasic
"""
    DiffusionSteadyDiph(phase1::Phase, phase2::Phase, bc_b::BorderConditions, ic::InterfaceConditions)

Creates a solver for a steady-state two-phase diffusion problem.

# Arguments
- `phase1::Phase`: The first phase of the problem.
- `phase2::Phase`: The second phase of the problem.
- `bc_b::BorderConditions`: The boundary conditions of the problem.
- `ic::InterfaceConditions`: The conditions at the interface between the two phases.
"""
function DiffusionSteadyDiph(phase1::Phase, phase2::Phase, bc_b::BorderConditions, ic::InterfaceConditions)
    println("Création du solveur:")
    println("- Diphasic problem")
    println("- Steady problem")
    println("- Diffusion problem")
    
    s = Solver(Steady, Diphasic, Diffusion, nothing, nothing, nothing, ConvergenceHistory(), [])
    
    s.A = build_diph_stead_diff_matrix(phase1.operator, phase2.operator, phase1.Diffusion_coeff, phase2.Diffusion_coeff, bc_b, ic)
    s.b = build_rhs_diph_stead_diff(phase1.operator, phase2.operator, phase1.source, phase2.source, phase1.capacity, phase2.capacity, bc_b, ic)

    BC_border_diph!(s.A, s.b, bc_b, phase2.capacity.mesh)

    return s
end

function build_diph_stead_diff_matrix(operator1::DiffusionOps, operator2::DiffusionOps, D1::Float64, D2::Float64, bc_b::BorderConditions, ic::InterfaceConditions)
    n = prod(operator1.size)

    jump, flux = ic.scalar, ic.flux
    Iₐ1, Iₐ2 = jump.α₁ * I(n), jump.α₂ * I(n)
    Iᵦ1, Iᵦ2 = flux.β₁ * I(n), flux.β₂ * I(n)
    Id1, Id2 = build_I_D(operator1, D1), build_I_D(operator2, D2)

    block1 = Id1 * operator1.G' * operator1.Wꜝ * operator1.G
    block2 = Id1 * operator1.G' * operator1.Wꜝ * operator1.H
    block3 = Id2 * operator2.G' * operator2.Wꜝ * operator2.G
    block4 = Id2 * operator2.G' * operator2.Wꜝ * operator2.H
    block5 = operator1.H' * operator1.Wꜝ * operator1.G
    block6 = operator1.H' * operator1.Wꜝ * operator1.H 
    block7 = operator2.H' * operator2.Wꜝ * operator2.G
    block8 = operator2.H' * operator2.Wꜝ * operator2.H

    A = spzeros(Float64, 4n, 4n)

    @inbounds begin
        # Top-left blocks
        A[1:n, 1:n] += block1
        A[1:n, n+1:2n] += block2

        # Middle blocks
        A[n+1:2n, n+1:2n] += Iₐ1
        A[n+1:2n, 3n+1:4n] -= Iₐ2

        # Bottom-left blocks
        A[2n+1:3n, 2n+1:3n] += block3
        A[2n+1:3n, 3n+1:4n] += block4

        # Bottom blocks with Iᵦ
        A[3n+1:4n, 1:n] += Iᵦ1 * block5
        A[3n+1:4n, n+1:2n] += Iᵦ1 * block6
        A[3n+1:4n, 2n+1:3n] += Iᵦ2 * block7
        A[3n+1:4n, 3n+1:4n] += Iᵦ2 * block8
    end

    return A
end

function build_rhs_diph_stead_diff(operator1::DiffusionOps, operator2::DiffusionOps, f1, f2, capacite1::Capacity, capacite2::Capacity, bc_b::BorderConditions, ic::InterfaceConditions)
    N = prod(operator1.size)
    b = zeros(4N)

    jump, flux = ic.scalar, ic.flux
    Iᵧ1, Iᵧ2 = build_I_g(operator1), build_I_g(operator2) #capacite1.Γ, capacite2.Γ #
    gᵧ, hᵧ = build_g_g(operator1, jump, capacite1), build_g_g(operator2, flux, capacite2)

    fₒ1 = build_source(operator1, f1, capacite1)
    fₒ2 = build_source(operator2, f2, capacite2)

    # Build the right-hand side
    b = vcat(operator1.V*fₒ1, gᵧ, operator2.V*fₒ2, Iᵧ2*hᵧ)

    return b
end

function solve_DiffusionSteadyDiph!(s::Solver, phase1::Phase, phase2::Phase; method::Function = gmres, kwargs...)
    if s.A === nothing
        error("Solver is not initialized. Call a solver constructor first.")
    end

    n = Int(size(s.A, 1) / 4)  # For diphasic problem, the system size is 4n

    if method == \
        # Remove zero rows and columns for direct solver
        s.A, s.b, rows_idx, cols_idx = remove_zero_rows_cols!(s.A, s.b)
        # Solve the reduced system
        x_reduced = s.A \ s.b
        # Reconstruct the full solution vector
        s.x = zeros(4n)
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


# Diffusion - Unsteady - Monophasic
"""
    DiffusionUnsteadyMono(phase::Phase, bc_b::BorderConditions, bc_i::AbstractBoundary, Δt::Float64, Tₑ::Float64, Tᵢ::Vector{Float64})

Constructs a solver for the unsteady monophasic diffusion problem.

# Arguments
- `phase::Phase`: The phase object representing the physical properties of the system.
- `bc_b::BorderConditions`: The border conditions object representing the boundary conditions at the outer border.
- `bc_i::AbstractBoundary`: The boundary conditions object representing the boundary conditions at the inner border.
- `Δt::Float64`: The time step size.
- `Tₑ::Float64`: The final time.
- `Tᵢ::Vector{Float64}`: The initial temperature distribution.
"""
function DiffusionUnsteadyMono(phase::Phase, bc_b::BorderConditions, bc_i::AbstractBoundary, Δt::Float64, Tₑ::Float64, Tᵢ::Vector{Float64})
    println("Création du solveur:")
    println("- Monophasic problem")
    println("- Unsteady problem")
    println("- Diffusion problem")
    
    s = Solver(Unsteady, Monophasic, Diffusion, nothing, nothing, nothing, ConvergenceHistory(), [])
    
    s.A = build_mono_unstead_diff_matrix(phase.operator, phase.capacity, phase.Diffusion_coeff, bc_b, bc_i, Δt)
    s.b = build_rhs_mono_unstead_diff(phase.operator, phase.source, phase.capacity, bc_b, bc_i, Tᵢ, Δt, 0.0)

    return s
end

function build_mono_unstead_diff_matrix(operator::DiffusionOps, capacite::Capacity, D::Float64, bc_b::BorderConditions, bc::AbstractBoundary, Δt::Float64)
    n = prod(operator.size)
    Iₐ, Iᵦ = build_I_bc(operator, bc)
    Iᵧ = build_I_g(operator) # capacite.Γ #
    Id = build_I_D(operator, D)

    # Preallocate the sparse matrix A with 2n rows and 2n columns
    A = spzeros(Float64, 2n, 2n)

    # Compute blocks
    block1 = operator.V + Δt / 2 * (Id * operator.G' * operator.Wꜝ * operator.G)
    block2 = Δt / 2 * (Id * operator.G' * operator.Wꜝ * operator.H)
    block3 = Iᵦ * operator.H' * operator.Wꜝ * operator.G
    block4 = Iᵦ * operator.H' * operator.Wꜝ * operator.H + (Iₐ * Iᵧ)

    A[1:n, 1:n] = block1
    A[1:n, n+1:2n] = block2
    A[n+1:2n, 1:n] = block3
    A[n+1:2n, n+1:2n] = block4
    
    return A
end

function build_rhs_mono_unstead_diff(operator::DiffusionOps, f, capacite::Capacity, bc_b::BorderConditions, bc::AbstractBoundary, Tᵢ, Δt::Float64, t::Float64)
    N = prod(operator.size)
    b = zeros(2N)

    Iᵧ = build_I_g(operator) #capacite.Γ #
    fₒn, fₒn1 = build_source(operator, f, t, capacite), build_source(operator, f, t+Δt, capacite)
    gᵧ = build_g_g(operator, bc, capacite)

    Tₒ, Tᵧ = Tᵢ[1:N], Tᵢ[N+1:end]

    # Build the right-hand side
    b = vcat((operator.V - Δt/2 * operator.G' * operator.Wꜝ * operator.G)*Tₒ - Δt/2 * operator.G' * operator.Wꜝ * operator.H * Tᵧ + Δt/2 * operator.V * (fₒn + fₒn1), Iᵧ * gᵧ)

   return b
end

function solve_DiffusionUnsteadyMono!(s::Solver, phase::Phase, Tᵢ, Δt::Float64, Tₑ, bc_b::BorderConditions, bc::AbstractBoundary; method::Function = gmres, kwargs...)
    if s.A === nothing
        error("Solver is not initialized. Call a solver constructor first.")
    end

    n = Int(size(s.A, 1) / 2)  # For monophasic problem, the system size is 2n

    t = 0.0
    while t < Tₑ
        t += Δt
        println("Time: ", t)
        s.b = build_rhs_mono_unstead_diff(phase.operator, phase.source, phase.capacity, bc_b, bc, Tᵢ, Δt, t)
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


# Diffusion - Unsteady - Diphasic
"""
    DiffusionUnsteadyDiph(phase1::Phase, phase2::Phase, bc_b::BorderConditions, ic::InterfaceConditions, Δt::Float64, Tₑ::Float64, Tᵢ::Vector{Float64})

Creates a solver for an unsteady two-phase diffusion problem.

## Arguments
- `phase1::Phase`: The first phase of the problem.
- `phase2::Phase`: The second phase of the problem.
- `bc_b::BorderConditions`: The boundary conditions of the problem.
- `ic::InterfaceConditions`: The conditions at the interface between the two phases.
- `Δt::Float64`: The time interval.
- `Tₑ::Float64`: The final time.
- `Tᵢ::Vector{Float64}`: The vector of initial temperatures.
"""
function DiffusionUnsteadyDiph(phase1::Phase, phase2::Phase, bc_b::BorderConditions, ic::InterfaceConditions, Δt::Float64, Tₑ::Float64, Tᵢ::Vector{Float64})
    println("Création du solveur:")
    println("- Diphasic problem")
    println("- Unsteady problem")
    println("- Diffusion problem")
    
    s = Solver(Unsteady, Diphasic, Diffusion, nothing, nothing, nothing, ConvergenceHistory(), [])
    
    s.A = build_diph_unstead_diff_matrix(phase1.operator, phase2.operator, phase1.Diffusion_coeff, phase2.Diffusion_coeff, bc_b, ic, Δt)
    s.b = build_rhs_diph_unstead_diff(phase1.operator, phase2.operator, phase1.source, phase2.source, phase1.capacity, phase2.capacity, bc_b, ic, Tᵢ, Δt, 0.0)

    return s
end

function build_diph_unstead_diff_matrix(operator1::DiffusionOps, operator2::DiffusionOps, D1::Float64, D2::Float64, bc_b::BorderConditions, ic::InterfaceConditions, Δt::Float64)
    n = prod(operator1.size)

    jump, flux = ic.scalar, ic.flux
    Iₐ1, Iₐ2 = jump.α₁ * I(n), jump.α₂ * I(n)
    Iᵦ1, Iᵦ2 = flux.β₁ * I(n), flux.β₂ * I(n)
    Id1, Id2 = build_I_D(operator1, D1), build_I_D(operator2, D2)

    # Precompute repeated multiplications
    WG_G1 = operator1.Wꜝ * operator1.G
    WG_H1 = operator1.Wꜝ * operator1.H
    WG_G2 = operator2.Wꜝ * operator2.G
    WG_H2 = operator2.Wꜝ * operator2.H

    block1 = operator1.V + Δt / 2 * Id1 * operator1.G' * WG_G1
    block2 = Δt / 2 * Id1 * operator1.G' * WG_H1
    block3 = operator2.V + Δt / 2 * Id2 * operator2.G' * WG_G2
    block4 = Δt / 2 * Id2 * operator2.G' * WG_H2
    block5 = Iᵦ1 * operator1.H' * WG_G1
    block6 = Iᵦ1 * operator1.H' * WG_H1
    block7 = Iᵦ2 * operator2.H' * WG_G2
    block8 = Iᵦ2 * operator2.H' * WG_H2

    # Preallocate the sparse matrix
    A = spzeros(Float64, 4n, 4n)

    # Assign blocks to the matrix
    A[1:n, 1:n] = block1
    A[1:n, n+1:2n] = block2
    A[1:n, 2n+1:3n] = spzeros(n, n)
    A[1:n, 3n+1:4n] = spzeros(n, n)

    A[n+1:2n, 1:n] = spzeros(n, n)
    A[n+1:2n, n+1:2n] = Iₐ1
    A[n+1:2n, 2n+1:3n] = spzeros(n, n)
    A[n+1:2n, 3n+1:4n] = -Iₐ2

    A[2n+1:3n, 1:n] = spzeros(n, n)
    A[2n+1:3n, n+1:2n] = spzeros(n, n)
    A[2n+1:3n, 2n+1:3n] = block3
    A[2n+1:3n, 3n+1:4n] = block4

    A[3n+1:4n, 1:n] = block5
    A[3n+1:4n, n+1:2n] = block6
    A[3n+1:4n, 2n+1:3n] = block7
    A[3n+1:4n, 3n+1:4n] = block8

    return A
end

function build_rhs_diph_unstead_diff(operator1::DiffusionOps, operator2::DiffusionOps, f1, f2, capacite1::Capacity, capacite2::Capacity, bc_b::BorderConditions, ic::InterfaceConditions, Tᵢ, Δt::Float64, t::Float64)
    N = prod(operator1.size)
    b = zeros(4N)

    jump, flux = ic.scalar, ic.flux
    Iᵧ1, Iᵧ2 = build_I_g(operator1), build_I_g(operator2) #capacite1.Γ, capacite2.Γ #
    gᵧ, hᵧ = build_g_g(operator1, jump,capacite1), build_g_g(operator2, flux, capacite2)

    fₒn1, fₒn2 = build_source(operator1, f1, t, capacite1), build_source(operator2, f2, t, capacite2)
    fₒn1p1, fₒn2p1 = build_source(operator1, f1, t+Δt, capacite1), build_source(operator2, f2, t+Δt, capacite2)

    Tₒ1, Tᵧ1 = Tᵢ[1:N], Tᵢ[N+1:2N]
    Tₒ2, Tᵧ2 = Tᵢ[2N+1:3N], Tᵢ[3N+1:end]

    # Build the right-hand side
    b = vcat((operator1.V - Δt/2 * operator1.G' * operator1.Wꜝ * operator1.G)*Tₒ1 - Δt/2 * operator1.G' * operator1.Wꜝ * operator1.H * Tᵧ1 + Δt/2 * operator1.V * (fₒn1 + fₒn1p1), gᵧ, (operator2.V - Δt/2 * operator2.G' * operator2.Wꜝ * operator2.G)*Tₒ2 - Δt/2 * operator2.G' * operator2.Wꜝ * operator2.H * Tᵧ2 + Δt/2 * operator2.V * (fₒn2 + fₒn2p1), Iᵧ2*hᵧ)

    return b
end

function solve_DiffusionUnsteadyDiph!(s::Solver, phase1::Phase, phase2::Phase, Tᵢ, Δt::Float64, Tₑ::Float64, bc_b::BorderConditions, ic::InterfaceConditions; method::Function = gmres, kwargs...)
    if s.A === nothing
        error("Solver is not initialized. Call a solver constructor first.")
    end

    n = Int(size(s.A, 1) / 4)  # For diphasic problem, the system size is 4n

    t = 0.0
    while t < Tₑ
        t += Δt
        println("Time: ", t)
        s.b = build_rhs_diph_unstead_diff(phase1.operator, phase2.operator, phase1.source, phase2.source, phase1.capacity, phase2.capacity, bc_b, ic, Tᵢ, Δt, t)
        BC_border_diph!(s.A, s.b, bc_b, phase2.capacity.mesh)
        
        if method == \
            # Remove zero rows and columns for direct solver
            A_reduced, b_reduced, _, cols_idx = remove_zero_rows_cols!(s.A, s.b)
            # Solve the reduced system
            x_reduced = A_reduced \ b_reduced
            # Reconstruct the full solution vector
            s.x = zeros(4n)
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
        Tᵢ = s.x
    end
end








# Advection - Unsteady - Monophasic
"""
    AdvectionUnsteadyMono(phase::Phase, bc_b::BorderConditions, bc_i::AbstractBoundary)

Creates a solver for a unsteady monophasic advection problem.

# Arguments
- `phase::Phase`: The phase object representing the physical properties of the system.
- `bc_b::BorderConditions`: The border conditions object representing the boundary conditions at the outer border.
- `bc_i::AbstractBoundary`: The boundary conditions object representing the boundary conditions at the inner border.
"""
function AdvectionUnsteadyMono(phase::Phase, bc_b::BorderConditions, bc_i::AbstractBoundary, Δt::Float64, Tₑ::Float64, Tᵢ::Vector{Float64})
    println("Création du solveur:")
    println("- Monophasic problem")
    println("- Unsteady problem")
    println("- Advection problem")
    
    s = Solver(Unsteady, Monophasic, Advection, nothing, nothing, nothing, ConvergenceHistory(), [])
    
    s.A = build_mono_unstead_adv_matrix(phase.operator, phase.capacity, bc_b, bc_i, Δt)
    s.b = build_rhs_mono_unstead_adv(phase.operator, phase.source, phase.capacity, bc_b, bc_i, Tᵢ, Δt, 0.0)

    return s
end

function build_mono_unstead_adv_matrix(operator::ConvectionOps, capacite::Capacity,  bc_b::BorderConditions, bc::AbstractBoundary, Δt::Float64)
    n = prod(operator.size)
    Iₐ, Iᵦ = build_I_bc(operator, bc)
    Iᵧ = build_I_g(operator) #capacite.Γ #

    C = operator.C # NTuple{N, SparseMatrixCSC{Float64, Int}}
    K = operator.K # NTuple{N, SparseMatrixCSC{Float64, Int}}

    A11 = operator.V + Δt/2 * (sum(C) + 0.5 * sum(K))
    A12 = +Δt/2 * 0.5 * sum(K)
    A21 = 0*Iᵦ # Iᵦ * operator.H' * operator.Wꜝ * operator.G
    A22 = 0*Iᵦ # Iᵦ * operator.H' * operator.Wꜝ * operator.H + Iₐ * Iᵧ

    A = vcat(hcat(A11, A12), hcat(A21, A22))
    return A
end

function build_rhs_mono_unstead_adv(operator::ConvectionOps, f, capacite::Capacity, bc_b::BorderConditions, bc::AbstractBoundary, Tᵢ, Δt::Float64, t::Float64)
    N = prod(operator.size)
    b = zeros(2N)

    C = operator.C # NTuple{N, SparseMatrixCSC{Float64, Int}}
    K = operator.K # NTuple{N, SparseMatrixCSC{Float64, Int}}

    Iᵧ = build_I_g(operator) #capacite.Γ #
    fₒn, fₒn1 = build_source(operator, f, t, capacite), build_source(operator, f, t+Δt, capacite)
    gᵧ = build_g_g(operator, bc, capacite)

    Tₒ, Tᵧ = Tᵢ[1:N], Tᵢ[N+1:end]

    # Build the right-hand side
    b1 = (operator.V - Δt/2 * sum(C) - Δt/2 * 0.5 * sum(K))*Tₒ - Δt/2 * 0.5 * sum(K) * Tᵧ + Δt/2 * operator.V * (fₒn + fₒn1)
    b2 = 0*Iᵧ * gᵧ
    b = vcat(b1, b2)
    return b
end

function solve_AdvectionUnsteadyMono!(s::Solver, phase::Phase, Tᵢ, Δt::Float64, Tₑ, bc_b::BorderConditions, bc::AbstractBoundary; method::Function = gmres, kwargs...)
    if s.A === nothing
        error("Solver is not initialized. Call a solver constructor first.")
    end

    n = Int(size(s.A, 1) / 2)  # For monophasic problem, the system size is 2n

    t = 0.0
    while t < Tₑ
        t += Δt
        println("Time: ", t)
        s.b = build_rhs_mono_unstead_adv(phase.operator, phase.source, phase.capacity, bc_b, bc, Tᵢ, Δt, t)
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








# AdvectionDiffusion - Steady - Monophasic
"""
    AdvectionDiffusionSteadyMono(phase::Phase, bc_b::BorderConditions, bc_i::AbstractBoundary)

Creates a solver for a steady-state monophasic advection-diffusion problem.

# Arguments
- `phase::Phase`: The phase object representing the physical properties of the system.
- `bc_b::BorderConditions`: The border conditions object representing the boundary conditions at the outer border.
- `bc_i::AbstractBoundary`: The boundary conditions object representing the boundary conditions at the inner border.
"""
function AdvectionDiffusionSteadyMono(phase::Phase, bc_b::BorderConditions, bc_i::AbstractBoundary)
    println("Création du solveur:")
    println("- Monophasic problem")
    println("- Steady problem")
    println("- Advection-Diffusion problem")
    
    s = Solver(Steady, Monophasic, DiffusionAdvection, nothing, nothing, nothing, ConvergenceHistory(), [])
    
    s.A = build_mono_stead_adv_diff_matrix(phase.operator, phase.capacity, phase.Diffusion_coeff, bc_b, bc_i)
    s.b = build_rhs_mono_stead_adv_diff(phase.operator, phase.source, phase.capacity, bc_b, bc_i)

    BC_border_mono!(s.A, s.b, bc_b, phase.capacity.mesh)

    return s
end

function build_mono_stead_adv_diff_matrix(operator::ConvectionOps, capacite::Capacity, D::Float64, bc_b::BorderConditions, bc::AbstractBoundary)
    n = prod(operator.size)
    Iₐ, Iᵦ = build_I_bc(operator, bc)
    Iᵧ = build_I_g(operator) #capacite.Γ #
    Id = build_I_D(operator, D)

    C = operator.C # NTuple{N, SparseMatrixCSC{Float64, Int}}
    K = operator.K # NTuple{N, SparseMatrixCSC{Float64, Int}}

    A11 = (sum(C) + 0.5 * sum(K)) + Id * operator.G' * operator.Wꜝ * operator.G
    A12 = 0.5 * sum(K) + Id * operator.G' * operator.Wꜝ * operator.H
    A21 = Iᵦ * operator.H' * operator.Wꜝ * operator.G
    A22 = Iᵦ * operator.H' * operator.Wꜝ * operator.H + Iₐ * Iᵧ

    A = vcat(hcat(A11, A12), hcat(A21, A22))
    return A
end

function build_rhs_mono_stead_adv_diff(operator::ConvectionOps, f, capacite::Capacity, bc_b::BorderConditions, bc::AbstractBoundary)
    N = prod(operator.size)
    b = zeros(2N)

    Iᵧ = build_I_g(operator) #capacite.Γ #
    fₒ = build_source(operator, f, capacite)
    gᵧ = build_g_g(operator, bc, capacite)

    # Build the right-hand side
    b = vcat(operator.V*fₒ, Iᵧ * gᵧ)

    return b
end

function solve_AdvectionDiffusionSteadyMono!(s::Solver, phase::Phase; method::Function = gmres, kwargs...)
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


# AdvectionDiffusion - Steady - Diphasic
"""
    AdvectionDiffusionSteadyDiph(phase1::Phase, phase2::Phase, bc_b::BorderConditions, ic::InterfaceConditions)

Creates a solver for a steady-state two-phase advection-diffusion problem.

# Arguments
- `phase1::Phase`: The first phase of the problem.
- `phase2::Phase`: The second phase of the problem.
- `bc_b::BorderConditions`: The boundary conditions of the problem.
- `ic::InterfaceConditions`: The conditions at the interface between the two phases.
"""
function AdvectionDiffusionSteadyDiph(phase1::Phase, phase2::Phase, bc_b::BorderConditions, ic::InterfaceConditions)
    println("Création du solveur:")
    println("- Diphasic problem")
    println("- Steady problem")
    println("- Advection-Diffusion problem")
    
    s = Solver(Steady, Diphasic, DiffusionAdvection, nothing, nothing, nothing, ConvergenceHistory(), [])
    
    s.A = build_diph_stead_adv_diff_matrix(phase1.operator, phase2.operator, phase1.Diffusion_coeff, phase2.Diffusion_coeff, bc_b, ic)
    s.b = build_rhs_diph_stead_adv_diff(phase1.operator, phase2.operator, phase1.source, phase2.source, phase1.capacity, phase2.capacity, bc_b, ic)

    BC_border_diph!(s.A, s.b, bc_b, phase2.capacity.mesh)

    return s
end

function build_diph_stead_adv_diff_matrix(operator1::ConvectionOps, operator2::ConvectionOps, D1::Float64, D2::Float64, bc_b::BorderConditions, ic::InterfaceConditions)
    n = prod(operator1.size)

    jump, flux = ic.scalar, ic.flux
    Iₐ1, Iₐ2 = jump.α₁*I(n), jump.α₂*I(n)
    Iᵦ1, Iᵦ2 = flux.β₁*I(n), flux.β₂*I(n)
    Id1, Id2 = build_I_D(operator1, D1), build_I_D(operator2, D2)

    C1 = operator1.C # NTuple{N, SparseMatrixCSC{Float64, Int}}
    K1 = operator1.K # NTuple{N, SparseMatrixCSC{Float64, Int}}
    C2 = operator2.C # NTuple{N, SparseMatrixCSC{Float64, Int}}
    K2 = operator2.K # NTuple{N, SparseMatrixCSC{Float64, Int}

    block1 = Id1 * operator1.G' * operator1.Wꜝ * operator1.G + (sum(C1) + 0.5 * sum(K1))
    block2 = Id1 * operator1.G' * operator1.Wꜝ * operator1.H + 0.5 * sum(K1)
    block3 = Id2 * operator2.G' * operator2.Wꜝ * operator2.G + (sum(C2) + 0.5 * sum(K2))
    block4 = Id2 * operator2.G' * operator2.Wꜝ * operator2.H + 0.5 * sum(K2)
    block5 = operator1.H' * operator1.Wꜝ * operator1.G
    block6 = operator1.H' * operator1.Wꜝ * operator1.H 
    block7 = operator2.H' * operator2.Wꜝ * operator2.G
    block8 = operator2.H' * operator2.Wꜝ * operator2.H

    A = vcat(hcat(block1, block2, zeros(n, n), zeros(n, n)),
             hcat(zeros(n, n), Iₐ1, zeros(n, n), -Iₐ2),
             hcat(zeros(n, n), zeros(n, n), block3, block4),
             hcat(Iᵦ1*block5, Iᵦ1*block6, Iᵦ2*block7, Iᵦ2*block8))
    return A
end

function build_rhs_diph_stead_adv_diff(operator1::ConvectionOps, operator2::ConvectionOps, f1, f2, capacite1::Capacity, capacite2::Capacity, bc_b::BorderConditions, ic::InterfaceConditions)
    N = prod(operator1.size)
    b = zeros(4N)

    jump, flux = ic.scalar, ic.flux
    Iᵧ1, Iᵧ2 = build_I_g(operator1), build_I_g(operator2) #capacite1.Γ, capacite2.Γ #
    gᵧ, hᵧ = build_g_g(operator1, jump, capacite1), build_g_g(operator2, flux, capacite2)

    fₒ1 = build_source(operator1, f1, capacite1)
    fₒ2 = build_source(operator2, f2, capacite2)

    # Build the right-hand side
    b = vcat(operator1.V*fₒ1, gᵧ, operator2.V*fₒ2, Iᵧ2*hᵧ)

    return b
end

function solve_AdvectionDiffusionSteadyDiph!(s::Solver, phase1::Phase, phase2::Phase; method::Function = gmres, kwargs...)
    if s.A === nothing
        error("Solver is not initialized. Call a solver constructor first.")
    end

    n = Int(size(s.A, 1) / 4)  # For diphasic problem, the system size is 4n

    if method == \
        # Remove zero rows and columns for direct solver
        s.A, s.b, rows_idx, cols_idx = remove_zero_rows_cols!(s.A, s.b)
        # Solve the reduced system
        x_reduced = s.A \ s.b
        # Reconstruct the full solution vector
        s.x = zeros(4n)
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


# AdvectionDiffusion - Unsteady - Monophasic
"""
    AdvectionDiffusionUnsteadyMono(phase::Phase, bc_b::BorderConditions, bc_i::AbstractBoundary, Δt::Float64, Tₑ::Float64, Tᵢ::Vector{Float64})

Creates a solver for an unsteady monophasic advection-diffusion problem.

# Arguments
- `phase::Phase`: The phase object representing the physical properties of the system.
- `bc_b::BorderConditions`: The border conditions object representing the boundary conditions at the outer border.
- `bc_i::AbstractBoundary`: The boundary conditions object representing the boundary conditions at the inner border.
- `Δt::Float64`: The time step size.
- `Tₑ::Float64`: The final time.
- `Tᵢ::Vector{Float64}`: The initial temperature distribution.
"""
function AdvectionDiffusionUnsteadyMono(phase::Phase, bc_b::BorderConditions, bc_i::AbstractBoundary, Δt::Float64, Tₑ::Float64, Tᵢ::Vector{Float64})
    println("Création du solveur:")
    println("- Monophasic problem")
    println("- Unsteady problem")
    println("- Advection-Diffusion problem")
    
    s = Solver(Unsteady, Monophasic, DiffusionAdvection, nothing, nothing, nothing, ConvergenceHistory(), [])
    
    s.A = build_mono_unstead_adv_diff_matrix(phase.operator, phase.capacity, phase.Diffusion_coeff, bc_b, bc_i, Δt)
    s.b = build_rhs_mono_unstead_adv_diff(phase.operator, phase.source, phase.capacity, bc_b, bc_i, Tᵢ, Δt, 0.0)

    return s
end

function build_mono_unstead_adv_diff_matrix(operator::ConvectionOps, capacite::Capacity, D::Float64, bc_b::BorderConditions, bc::AbstractBoundary, Δt::Float64)
    n = prod(operator.size)
    Iₐ, Iᵦ = build_I_bc(operator, bc)
    Iᵧ = build_I_g(operator) #capacite.Γ #
    Id = build_I_D(operator, D)

    C = operator.C # NTuple{N, SparseMatrixCSC{Float64, Int}}
    K = operator.K # NTuple{N, SparseMatrixCSC{Float64, Int}}

    A11 = operator.V + Δt/2 * (sum(C) + 0.5 * sum(K)) + Δt/2 * Id * operator.G' * operator.Wꜝ * operator.G
    A12 = Δt/2 * 0.5 * sum(K) + Δt/2 * Id * operator.G' * operator.Wꜝ * operator.H
    A21 = Iᵦ * operator.H' * operator.Wꜝ * operator.G
    A22 = Iᵦ * operator.H' * operator.Wꜝ * operator.H + Iₐ * Iᵧ

    A = vcat(hcat(A11, A12), hcat(A21, A22))
    return A
end

function build_rhs_mono_unstead_adv_diff(operator::ConvectionOps, f, capacite::Capacity, bc_b::BorderConditions, bc::AbstractBoundary, Tᵢ, Δt::Float64, t::Float64)
    N = prod(operator.size)
    b = zeros(2N)

    C = operator.C # NTuple{N, SparseMatrixCSC{Float64, Int}}
    K = operator.K # NTuple{N, SparseMatrixCSC{Float64, Int}}

    Iᵧ = build_I_g(operator) #capacite.Γ #
    fₒn, fₒn1 = build_source(operator, f, t, capacite), build_source(operator, f, t+Δt, capacite)
    gᵧ = build_g_g(operator, bc, capacite)

    Tₒ, Tᵧ = Tᵢ[1:N], Tᵢ[N+1:end]

    """
    lx, ly = 16., 16.
    nx, ny = 160, 160
    radius, center = ly/4, (lx/2, ly/2) .+ (0.01, 0.01)

    # Coordinates of the grid centers
    x_coords = capacite.mesh.centers[1]
    y_coords = capacite.mesh.centers[2]

    # Define the radius of the circle (in physical units)
    circle_radius = lx / 30  # Adjust this value to control the size of the circle

    # Loop over all grid points
    for j in 1:(ny )
        for i in 1:(nx )
            idx = i + (j - 1) * (nx + 1)
            x_i = x_coords[i]
            y_j = y_coords[j]
            # Compute the distance from the center
            distance = sqrt((x_i - center[1])^2 + (y_j - center[2])^2)
            # If the point is inside the circle, set T=1
            if distance <= circle_radius
                Tₒ[idx] = 1.0
                Tᵧ[idx] = 1.0
            end
        end
    end
    """

    # Build the right-hand side
    b = vcat(operator.V * Tₒ  - Δt/2 * sum(C) * Tₒ - 0.5 * sum(K) * Tₒ - Δt/2 * 0.5 * sum(K) * Tᵧ - Δt/2 * operator.G' * operator.Wꜝ * operator.G * Tₒ - Δt/2 * operator.G' * operator.Wꜝ * operator.H * Tᵧ + Δt/2 * operator.V * (fₒn + fₒn1), Iᵧ * gᵧ)
    return b
end

function solve_AdvectionDiffusionUnsteadyMono!(s::Solver, phase::Phase, Tᵢ, Δt::Float64, Tₑ, bc_b::BorderConditions, bc::AbstractBoundary; method::Function = gmres, kwargs...)
    if s.A === nothing
        error("Solver is not initialized. Call a solver constructor first.")
    end

    n = Int(size(s.A, 1) / 2)  # For monophasic problem, the system size is 2n

    t = 0.0
    while t < Tₑ
        t += Δt
        println("Time: ", t)
        s.b = build_rhs_mono_unstead_adv_diff(phase.operator, phase.source, phase.capacity, bc_b, bc, Tᵢ, Δt, t)
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


# AdvectionDiffusion - Unsteady - Diphasic
"""
    AdvectionDiffusionUnsteadyDiph(phase1::Phase, phase2::Phase, bc_b::BorderConditions, ic::InterfaceConditions, Δt::Float64, Tₑ::Float64, Tᵢ::Vector{Float64})

Creates a solver for an unsteady two-phase advection-diffusion problem.

# Arguments
- `phase1::Phase`: The first phase of the problem.
- `phase2::Phase`: The second phase of the problem.
- `bc_b::BorderConditions`: The boundary conditions of the problem.
- `ic::InterfaceConditions`: The conditions at the interface between the two phases.
- `Δt::Float64`: The time interval.
- `Tₑ::Float64`: The final time.
- `Tᵢ::Vector{Float64}`: The vector of initial temperatures.
"""
function AdvectionDiffusionUnsteadyDiph(phase1::Phase, phase2::Phase, bc_b::BorderConditions, ic::InterfaceConditions, Δt::Float64, Tₑ::Float64, Tᵢ::Vector{Float64})
    println("Création du solveur:")
    println("- Diphasic problem")
    println("- Unsteady problem")
    println("- Advection-Diffusion problem")
    
    s = Solver(Unsteady, Diphasic, DiffusionAdvection, nothing, nothing, nothing, ConvergenceHistory(), [])
    
    s.A = build_diph_unstead_adv_diff_matrix(phase1.operator, phase2.operator, phase1.Diffusion_coeff, phase2.Diffusion_coeff, bc_b, ic, Δt)
    s.b = build_rhs_diph_unstead_adv_diff(phase1.operator, phase2.operator, phase1.source, phase2.source, phase1.capacity, phase2.capacity, bc_b, ic, Tᵢ, Δt, 0.0)

    return s
end

function build_diph_unstead_adv_diff_matrix(operator1::ConvectionOps, operator2::ConvectionOps, D1::Float64, D2::Float64, bc_b::BorderConditions, ic::InterfaceConditions, Δt::Float64)
    n = prod(operator1.size)

    jump, flux = ic.scalar, ic.flux
    Iₐ1, Iₐ2 = jump.α₁*I(n), jump.α₂*I(n)
    Iᵦ1, Iᵦ2 = flux.β₁*I(n), flux.β₂*I(n)
    Id1, Id2 = build_I_D(operator1, D1), build_I_D(operator2, D2)

    C1 = operator1.C # NTuple{N, SparseMatrixCSC{Float64, Int}}
    K1 = operator1.K # NTuple{N, SparseMatrixCSC{Float64, Int}}
    C2 = operator2.C # NTuple{N, SparseMatrixCSC{Float64, Int}}
    K2 = operator2.K # NTuple{N, SparseMatrixCSC{Float64, Int}

    block1 = operator1.V + Δt/2 * (sum(C1) + 0.5 * sum(K1)) + Δt/2 * Id1 * operator1.G' * operator1.Wꜝ * operator1.G
    block2 = Δt/2 * 0.5 * sum(K1) + Δt/2 * Id1 * operator1.G' * operator1.Wꜝ * operator1.H
    block3 = operator2.V + Δt/2 * (sum(C2) + 0.5 * sum(K2)) + Δt/2 * Id2 * operator2.G' * operator2.Wꜝ * operator2.G
    block4 = Δt/2 * 0.5 * sum(K2) + Δt/2 * Id2 * operator2.G' * operator2.Wꜝ * operator2.H
    block5 = Iᵦ1 * operator1.H' * operator1.Wꜝ * operator1.G
    block6 = Iᵦ1 * operator1.H' * operator1.Wꜝ * operator1.H
    block7 = Iᵦ2 * operator2.H' * operator2.Wꜝ * operator2.G
    block8 = Iᵦ2 * operator2.H' * operator2.Wꜝ * operator2.H

    A = vcat(hcat(block1, block2, zeros(n, n), zeros(n, n)),
             hcat(zeros(n, n), Iₐ1, zeros(n, n), -Iₐ2),
             hcat(zeros(n, n), zeros(n, n), block3, block4),
             hcat(Iᵦ1*block5, Iᵦ1*block6, Iᵦ2*block7, Iᵦ2*block8))

    return A
end

function build_rhs_diph_unstead_adv_diff(operator1::ConvectionOps, operator2::ConvectionOps, f1, f2, capacite1::Capacity, capacite2::Capacity, bc_b::BorderConditions, ic::InterfaceConditions, Tᵢ, Δt::Float64, t::Float64)
    N = prod(operator1.size)
    b = zeros(4N)

    jump, flux = ic.scalar, ic.flux
    Iᵧ1, Iᵧ2 = build_I_g(operator1), build_I_g(operator2) #capacite1.Γ, capacite2.Γ #
    gᵧ, hᵧ = build_g_g(operator1, jump,capacite1), build_g_g(operator2, flux, capacite2)

    fₒn1, fₒn2 = build_source(operator1, f1, t, capacite1), build_source(operator2, f2, t, capacite2)
    fₒn1p1, fₒn2p1 = build_source(operator1, f1, t+Δt, capacite1), build_source(operator2, f2, t+Δt, capacite2)

    Tₒ1, Tᵧ1 = Tᵢ[1:N], Tᵢ[N+1:2N]
    Tₒ2, Tᵧ2 = Tᵢ[2N+1:3N], Tᵢ[3N+1:end]

    # Build the right-hand side
    b = vcat((operator1.V - Δt/2 * sum(operator1.C) - Δt/2 * 0.5 * sum(operator1.K))*Tₒ1 - Δt/2 * 0.5 * sum(operator1.K) * Tᵧ1 + Δt/2 * operator1.V * (fₒn1 + fₒn1p1), gᵧ, (operator2.V - Δt/2 * sum(operator2.C) - Δt/2 * 0.5 * sum(operator2.K))*Tₒ2 - Δt/2 * 0.5 * sum(operator2.K) * Tᵧ2 + Δt/2 * operator2.V * (fₒn2 + fₒn2p1), Iᵧ2*hᵧ)

    return b
end

function solve_AdvectionDiffusionUnsteadyDiph!(s::Solver, phase1::Phase, phase2::Phase, Tᵢ, Δt::Float64, Tₑ, bc_b::BorderConditions, ic::InterfaceConditions; method::Function = gmres, kwargs...)
    if s.A === nothing
        error("Solver is not initialized. Call a solver constructor first.")
    end

    n = Int(size(s.A, 1) / 4)  # For diphasic problem, the system size is 4n

    t = 0.0
    while t < Tₑ
        t += Δt
        println("Time: ", t)
        s.b = build_rhs_diph_unstead_adv_diff(phase1.operator, phase2.operator, phase1.source, phase2.source, phase1.capacity, phase2.capacity, bc_b, ic, Tᵢ, Δt, t)
        BC_border_diph!(s.A, s.b, bc_b, phase2.capacity.mesh)
        
        if method == \
            # Remove zero rows and columns for direct solver
            A_reduced, b_reduced, _, cols_idx = remove_zero_rows_cols!(s.A, s.b)
            # Solve the reduced system
            x_reduced = A_reduced \ b_reduced
            # Reconstruct the full solution vector
            s.x = zeros(4n)
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
        Tᵢ = s.x
    end
end






function BC_border_mono!(A::SparseMatrixCSC{Float64, Int}, b::Vector{Float64}, bc_b::BorderConditions, mesh::CartesianMesh) 
    # Identify location border cells based on the dimension : 1D (left, right), 2D (left, right, top, bottom), 3D (left, right, top, bottom, front, back)
    left_cells = []
    right_cells = []
    top_cells = []
    bottom_cells = []
    forward_cells = []
    backward_cells = []
    for (key, bc) in bc_b.borders
        if key == :left
            left_cells = [cell for (i, (cell, linear_index)) in enumerate(mesh.tag.border_cells) if cell[2] == 1]
        elseif key == :right
            right_cells = [cell for (i, (cell, linear_index)) in enumerate(mesh.tag.border_cells) if cell[2] == length(mesh.centers[2])]
        elseif key == :top
            top_cells = [cell for (i, (cell, linear_index)) in enumerate(mesh.tag.border_cells) if cell[1] == length(mesh.centers[1])]
        elseif key == :bottom
            bottom_cells = [cell for (i, (cell, linear_index)) in enumerate(mesh.tag.border_cells) if cell[1] == 1]
        elseif key == :forward
            forward_cells = [cell for (i, (cell, linear_index)) in enumerate(mesh.tag.border_cells) if cell[3] == length(mesh.centers[3])]
        elseif key == :backward
            backward_cells = [cell for (i, (cell, linear_index)) in enumerate(mesh.tag.border_cells) if cell[3] == 1]
        end
    end
    for (i, (cell, linear_index)) in enumerate(mesh.tag.border_cells)
        condition = nothing
        current_key = nothing
        if cell in left_cells
            condition = bc_b.borders[:left]
            current_key = :left
        elseif cell in right_cells
            condition = bc_b.borders[:right]
            current_key = :right
        elseif cell in top_cells
            condition = bc_b.borders[:top]
            current_key = :top
        elseif cell in bottom_cells
            condition = bc_b.borders[:bottom]
            current_key = :bottom
        elseif cell in forward_cells
            condition = bc_b.borders[:forward]
            current_key = :forward
        elseif cell in backward_cells
            condition = bc_b.borders[:backward]
            current_key = :backward
        end

        if condition isa Dirichlet
            A[linear_index, :] .= 0.0
            A[linear_index, linear_index] = 1.0
            b[linear_index] = isa(condition.value, Function) ? condition.value(cell...) : condition.value
        elseif condition isa Periodic
            # Apply Periodic condition
            opposite_key = get_opposite_boundary(current_key)
            if !haskey(bc_b.borders, opposite_key)
                error("Periodic boundary requires both boundaries to be specified")
            end

            # Find corresponding cell index
            corresponding_cell = find_corresponding_cell(cell, current_key, opposite_key, mesh)
            corresponding_idx = cell_to_index(mesh, corresponding_cell)

            # Modify A to enforce x_i - x_j = 0
            A[linear_index, linear_index] += 1.0
            A[linear_index, corresponding_idx] -= 1.0
            b[linear_index] = 0.0
        elseif condition isa Neumann
            # Not implemented yet
        elseif condition isa Robin
            # Not implemented yet
        end
    end
end

# Helper function to get the opposite boundary
function get_opposite_boundary(key::Symbol)
    if key == :left
        return :right
    elseif key == :right
        return :left
    elseif key == :bottom
        return :top
    elseif key == :top
        return :bottom
    elseif key == :backward
        return :forward
    elseif key == :forward
        return :backward
    else
        error("Unknown boundary key: $key")
    end
end

# Helper function to find the corresponding cell on the opposite boundary
function find_corresponding_cell(cell::CartesianIndex{N}, key::Symbol, opposite_key::Symbol, mesh::CartesianMesh) where {N}
    if key == :left || key == :right
        new_cell = CartesianIndex(key == :left ? length(mesh.centers[1]) : 1, cell[2])
    elseif key == :bottom || key == :top
        new_cell = CartesianIndex(cell[1], key == :bottom ? length(mesh.centers[2]) : 1)
    elseif key == :backward || key == :forward
        new_cell = CartesianIndex(cell[1], cell[2], key == :backward ? length(mesh.centers[3]) : 1)
    end
    return new_cell
end

function cell_to_index(mesh::CartesianMesh{1}, cell::CartesianIndex)
    return LinearIndices((length(mesh.centers[1])+1,))[cell]
end

function cell_to_index(mesh::CartesianMesh{2}, cell::CartesianIndex)
    return LinearIndices((length(mesh.centers[1])+1, length(mesh.centers[2])+1))[cell]
end

function cell_to_index(mesh::CartesianMesh{3}, cell::CartesianIndex)
    return LinearIndices((length(mesh.centers[1])+1, length(mesh.centers[2])+1, length(mesh.centers[3])+1))[cell]
end

function BC_border_diph!(A::SparseMatrixCSC{Float64, Int}, b::Vector{Float64}, bc_b::BorderConditions, mesh::CartesianMesh) 
    # Identify location border cells based on the dimension : 1D (left, right), 2D (left, right, top, bottom), 3D (left, right, top, bottom, front, back)
    left_cells = []
    right_cells = []
    top_cells = []
    bottom_cells = []
    forward_cells = []
    backward_cells = []
    for (key, bc) in bc_b.borders
        if key == :left
            left_cells = [cell for (i, (cell, linear_index)) in enumerate(mesh.tag.border_cells) if cell[2] == 1]
        elseif key == :right
            right_cells = [cell for (i, (cell, linear_index)) in enumerate(mesh.tag.border_cells) if cell[2] == length(mesh.centers[2])]
        elseif key == :top
            top_cells = [cell for (i, (cell, linear_index)) in enumerate(mesh.tag.border_cells) if cell[1] == length(mesh.centers[1])]
        elseif key == :bottom
            bottom_cells = [cell for (i, (cell, linear_index)) in enumerate(mesh.tag.border_cells) if cell[1] == 1]
        elseif key == :forward
            forward_cells = [cell for (i, (cell, linear_index)) in enumerate(mesh.tag.border_cells) if cell[3] == length(mesh.centers[3])]
        elseif key == :backward
            backward_cells = [cell for (i, (cell, linear_index)) in enumerate(mesh.tag.border_cells) if cell[3] == 1]
        end
    end
    for (i, (cell, linear_index)) in enumerate(mesh.tag.border_cells)
        condition = nothing
        if cell in left_cells
            condition = bc_b.borders[:left]
        elseif cell in right_cells
            condition = bc_b.borders[:right]
        elseif cell in top_cells
            condition = bc_b.borders[:top]
        elseif cell in bottom_cells
            condition = bc_b.borders[:bottom]
        elseif cell in forward_cells
            condition = bc_b.borders[:forward]
        elseif cell in backward_cells
            condition = bc_b.borders[:backward]
        end

        if linear_index!=1
            linear_index = linear_index + size(A, 1) ÷ 2 
        end

        if condition isa Dirichlet
            A[linear_index, :] .= 0.0
            A[linear_index, linear_index] = 1.0
            b[linear_index] = isa(condition.value, Function) ? condition.value(cell...) : condition.value
        elseif condition isa Neumann
            # Not implemented yet
        elseif condition isa Robin
            # Not implemented yet
        end
    end
end
    
function build_I_D(operator::AbstractOperators, D::Union{Float64,Function})
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

function build_g_g(operator::AbstractOperators, bc::Union{AbstractBoundary, AbstractInterfaceBC}, capacite::Capacity)
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

function build_source(operator::AbstractOperators, f, capacite::Capacity)
    N = prod(operator.size)
    fₒ = zeros(N)

    # Compute the source term
    for i in 1:N
        x, y, z = get_coordinates(i, capacite.C_ω)
        fₒ[i] = f(x, y, z)
    end

    return fₒ
end

function build_source(operator::AbstractOperators, f, t, capacite::Capacity)
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
