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


# TO ADAPT
function measure(sdf, map, x, t)
    # Define a function that depends only on x, by setting t
    f(x_vec) = sdf(x_vec, t)
    
    # Calculate the signed distance
    d = f(x)
    
    # Calculate the gradient of the signed distance with respect to x
    ∇d = ForwardDiff.gradient(f, x)
    
    # Check if the gradient is zero to avoid division by zero
    if norm(∇d) == 0
        n = zeros(Float64, length(x))
    else
        # Calculate the normal vector
        n = ∇d ./ norm(∇d)
    end
    
    # Calculate the Jacobian matrix of the map function with respect to x
    J = ForwardDiff.jacobian(x_vec -> map(x_vec, t), x)
    
    # Calculate the time derivative of the map function
    dot = ForwardDiff.derivative(t_val -> map(x, t_val), t)
    
    # Calculate velocity
    v = -J \ dot
    
    return d, n, v
end
