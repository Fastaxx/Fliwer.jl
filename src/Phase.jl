struct Phase
    capacity::AbstractCapacity
    operator::AbstractOperators
    source::Union{Function, Nothing}
    Diffusion_coeff::Union{Function, Float64}
end


