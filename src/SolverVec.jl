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
function NavierStokesUnsteadyMono(pressure, velocity, bc_b, Δt, Tend, uₒ)
    println("Création du solveur Navier-Stokes Monophasic Unsteady")
    println("- Monophasic problem")
    println("- Unsteady problem")

    s = VectorSolver(Unsteady, Monophasic, DiffusionAdvection, nothing, nothing, nothing, ConvergenceHistory(), [])

    s.A = build_navier_stokes_matrix(pressure.capacity, velocity.capacity_u, velocity.capacity_v, velocity.ρ, velocity.Re, Δt)
    s.b = build_navier_stokes_rhs(pressure.capacity, velocity.capacity_u, velocity.capacity_v, velocity.ρ, velocity.Re, Δt, uₒ)

    # BC_border_mono!(s.A, s.b, bc_b, phase.capacity.mesh)
    return s
end

function build_navier_stokes_matrix(capacity_p, capacity_u, capacity_v, ρ, Re, Δt)
    return nothing
end

function build_navier_stokes_rhs(capacity_p, capacity_u, capacity_v, ρ, Re, Δt, uₒ)
    return nothing
end

function solve_NavierStokesUnsteadyMono!(solver::VectorSolver, pressure::Pressure, velocity::Velocity, uₒ, Δt, Tend, bc; method=IterativeSolvers.bicgstabl, abstol=1e-15, verbose=false)
    println("Solve Navier-Stokes Monophasic Unsteady")
    println("- Monophasic problem")
    println("- Unsteady problem")
    println("- Solver: ", method)
    println("- Absolute tolerance: ", abstol)
    println("- Verbose: ", verbose)
end
    