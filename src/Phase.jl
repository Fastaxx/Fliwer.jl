"""
    struct Phase

The `Phase` struct represents a phase in a system.

# Fields
- `capacity::AbstractCapacity`: The capacity of the phase.
- `operator::AbstractOperators`: The operators associated with the phase.
- `source::Union{Function, Nothing}`: The source function or `Nothing` if there is no source.
- `Diffusion_coeff::Union{Function, Float64}`: The diffusion coefficient function or a constant value.
"""
mutable struct Phase
    capacity::AbstractCapacity
    operator::AbstractOperators
    source::Union{Function, Nothing}
    Diffusion_coeff::Union{Function, Float64}
end

mutable struct Velocity
    capacity_u::AbstractCapacity
    capacity_v::AbstractCapacity
    operator::AbstractOperators
    source::Union{Function, Nothing}
    ρ::Float64
    Re::Float64
end

mutable struct Pressure
    capacity::AbstractCapacity
    operator::AbstractOperators
    source::Union{Function, Nothing}
    ρ::Float64
end

