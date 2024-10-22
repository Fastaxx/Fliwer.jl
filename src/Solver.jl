abstract type AbstractSolver end

struct DiffusionSolver{N} <: AbstractSolver where N
    ops::DiffusionOps{N}
    steady_state::Bool
    dt::Float64
    t::Float64
end
