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

"""
    d = sdf(body::AutoBody,x,t) = body.sdf(x,t)
"""
sdf(body::Body, x, t=0; kwargs...) = body.sdf(x, t; kwargs...)

Base.:+(a::Body, b::Body) = Body((x,t)->min(a.sdf(x,t), b.sdf(x,t)), (x,t)->ifelse(a.sdf(x,t) < b.sdf(x,t), a.map(x,t), b.map(x,t)), a.domain, a.compose || b.compose)
∩(a::Body, b::Body) = Body((x,t)->max(a.sdf(x,t), b.sdf(x,t)), (x,t)->ifelse(a.sdf(x,t) > b.sdf(x,t), a.map(x,t), b.map(x,t)), a.domain, a.compose || b.compose)
⊖(a::Body, b::Body) = Body((x,t)->max(a.sdf(x,t), -b.sdf(x,t)), a.map, a.domain, a.compose)
c(a::Body) = Body((x,t)->-a.sdf(x,t), a.map, a.domain, false)

measure(body::Body, x, t) = measure(body.sdf, body.map, x, t)

"""
    d, n, v = measure(sdf, map, x, t)

Determine the signed distance, normal vector, and velocity of a body at a given point.
"""
function measure(sdf, map, x, t)
    # Eval the sdf
    d = sdf(x...)

    # Eval the gradient of the sdf
    f = x -> sdf(x...)
    ∇d = ForwardDiff.gradient(f, [x...])

    # Check if the gradient is zero to avoid division by zero
    if norm(∇d) == 0
        n = zeros(Float64, length(x))
    else
        # Calculate the normal vector
        n = ∇d ./ norm(∇d)
    end

    # Calculate the Jacobian matrix of the map function with respect to x 
    # Dm/Dt=0 → ṁ + (dm/dx)ẋ = 0 ∴  ẋ =-(dm/dx)\ṁ
    J = ForwardDiff.jacobian(u -> [map(u..., t)[i] for i in 1:length(x)], x)
    dt = ForwardDiff.derivative(u -> sum([map(x..., u)...]), t) # A VERIFIER

    # Calculate velocity
    v = dt ./ J

    return d, n, v
end
