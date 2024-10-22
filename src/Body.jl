"""
    AbstractBody

Embedded body Abstract type. Any `AbstractBody` subtype must implement

d = sdf(body::AbstractBody, x, t=0)
"""
abstract type AbstractBody end

"""
    NoBody

Use for a simulation without any body.
"""
struct NoBody <: AbstractBody end

"""
    struct Body{N} <: AbstractBody

The `Body` struct represents a body in the Fliwer package.

# Fields
- `sdf::Function`: The signed distance function of the body.
- `map::Function`: The mapping function of the body.
- `domain::NTuple{N, Tuple{Float64, Float64}}`: The domain of the body.
- `compose::Bool`: Flag for composing `sdf=sdf∘map`

"""
struct Body{N} <: AbstractBody
    sdf::Function
    map::Function
    domain::NTuple{N, Tuple{Float64, Float64}}
    compose::Bool
end

sdf(body::Body, x, t=0; kwargs...) = body.sdf(x, t; kwargs...)
