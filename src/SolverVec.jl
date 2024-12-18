mutable struct VectorSolver{TT<:TimeType, PT<:PhaseType, ET<:EquationType}
    time_type::TT
    phase_type::PT
    equation_type::ET
    A::Union{SparseMatrixCSC{Float64, Int}, Nothing}
    b::Union{Vector{Float64}, Nothing}
    x::Union{Vector{Float64}, Nothing}
    ch::IterativeSolvers.ConvergenceHistory
    states::Vector{Any}
end

# Navier-Stokes Monophasic Unsteady
function NavierStokesUnsteadyMono(velocity, bc_b, Δt, Tend, x0)
    println("Création du solveur Navier-Stokes Monophasic Unsteady")
    println("- Monophasic problem")
    println("- Unsteady problem")

    s = VectorSolver(Unsteady, Monophasic, DiffusionAdvection, nothing, nothing, nothing, ConvergenceHistory(), [])

    if typeof(velocity)==Velocity{1}
        s.A = build_navier_stokes_matrix(velocity.operator, velocity.capacities, velocity.ρ, velocity.Re, Δt, x0)
        s.b = build_navier_stokes_rhs(velocity.operator, velocity.capacities, velocity.ρ, velocity.Re, Δt, x0)
    elseif typeof(velocity)==Velocity{2}
        println("2D")
    end

    # BC_border_mono!(s.A, s.b, bc_b, phase.capacity.mesh)
    return s
end

function build_navier_stokes_matrix(operator, capacity_u, ρ, Re, Δt, x0)
    n = prod(operator.size)

    #A = spzeros(Float64, n, n)
    Vx = operator.V[1]
    Gx = operator.G[1]
    Wꜝx= operator.Wꜝ[1]
    Hx = operator.H[1]

    # Compute the blocks
    block1 = Vx/Δt + 1/(2*Re) * Gx' * Wꜝx * Gx  # Add the convective term
    block2 = 1/(2*Re) * Gx' * Wꜝx * Hx # Add the convective term
    block3 = 1/ρ * Wꜝx * Gx
    block4 = 1/ρ * Wꜝx * Hx
    block5 = I(n)
    block6 = -(Gx' + Hx')
    block7 = Hx'

    # Fill the matrix
    A = [block1 block2 block3 block4;
         zeros(n, n) block5 zeros(n, n) zeros(n, n);
         block6 block7 zeros(n, n) zeros(n, n);
         zeros(n, n) zeros(n, n) zeros(n, n) zeros(n, n)]

    return A
end

function build_navier_stokes_rhs(operator, capacity_u, ρ, Re, Δt, x0)
    n = prod(operator.size)
    b = zeros(4n)

    Vx = operator.V[1]
    Gx = operator.G[1]
    Wꜝx= operator.Wꜝ[1]
    Hx = operator.H[1]
    Cx = operator.C[1]
    Kx = operator.K[1]

    uₒ = (x0[1:n],)
    uᵧ = x0[n+1:2n]

    # Compute the convective terms
    convxn = Cx * uₒ[1] + Kx * (uₒ[1] + uᵧ)/2  # Recompute with the new uₒ and uᵧ and the old (n-1)
    convxn1 = Cx * uₒ[1] + Kx * (uₒ[1] + uᵧ)/2

    b1 = Vx/Δt * uₒ[1] - 1/(2*Re) * Gx' * Wꜝx * Gx * uₒ[1] - 1/(2*Re) * Gx' * Wꜝx * Hx * uᵧ - 3/2 * convxn + 1/2 * convxn1
    b2 = zeros(n)
    b3 = zeros(n)
    b4 = zeros(n)

    b = [b1; b2; b3; b4]
    return b
end

function solve_NavierStokesUnsteadyMono!(s::VectorSolver, velocity::Velocity, uₒ, uᵧ, Δt, Tend, bc; method=IterativeSolvers.bicgstabl, kwargs...)
    if s.A === nothing
        error("Solver is not initialized. Call a solver constructor first.")
    end

    n = Int(size(s.A, 1) / 4) # For 1D

    t=0.0

    s.x = method(s.A, s.b; kwargs...)
    push!(s.states, s.x)
    uₒ = (s.x[1:n],)
    uᵧ = s.x[n+1:2n]
    while t<Tend
        t+=Δt
        println("t = ", t)
        # Update the right-hand side
        if typeof(velocity)==Velocity{1}
            s.b = build_navier_stokes_rhs(velocity.operator, velocity.capacities, velocity.ρ, velocity.Re, Δt, s.x)
        elseif typeof(velocity)==Velocity{2}
            println("2D")
        end

        #Apply the boundary conditions
        #BC_border_mono!(s.A, s.b, bc, phase.capacity.mesh)
        
        # Solve the linear system
        s.x = method(s.A, s.b; kwargs...)

        # Update the solution
        push!(s.states, s.x)
        uₒ = (s.x[1:n],)
        uᵧ = s.x[n+1:2n]
    end
end
    