@enum TimeType begin
    Steady
    Unsteady
end

@enum PhaseType begin
    Monophasic
    Diphasic
end

@enum EquationType begin
    Advection            # ∂ₜu⃗ +(u⃗ ·∇)u⃗ = 0
end

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

# Advection - Steady - Mono
function VectorAdvectionSteadyMono(phase::Phase, bc_b::BorderConditions, bc_i::AbstractBoundary)
    println("Création du solveur Vector")
    println("- Monophasic problem")
    println("- Steady problem")
    println("- Advection problem")

    s = VectorSolver(Steady, Monophasic, Advection, nothing, nothing, nothing, ConvergenceHistory(), [])

    s.A = build_vector_mono_steady_advection_matrix(phase.operator, phase.capacity, phase.Diffusion_coeff, bc_b, bc_i)
    s.b = build_vector_mono_steady_advection_rhs(phase.operator, phase.capacity, bc_b, bc_i)

    # BC_border_mono!(s.A, s.b, bc_b, phase.capacity.mesh)
    return s
end

function build_vector_mono_steady_advection_matrix(operator::AbstractOperators, capacity::Capacity, Diffusion_coeff::Float64, bc_b::BorderConditions, bc_i::AbstractBoundary)
    
end