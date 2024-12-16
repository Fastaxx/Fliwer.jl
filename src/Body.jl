"""
    AbstractBody

Embedded body Abstract type. Any `AbstractBody` subtype must implement

d = sdf(body::AbstractBody, x, t=0)
"""
abstract type AbstractBody end

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

"""
    NoBody

Use for a simulation without any body.
"""
struct NoBody <: AbstractBody end

NoBody1D(domain) = Body((x,_=0)->-1.0, (x,_)->(x,), domain, false)
NoBody2D(domain) = Body((x,y,_=0)->-1.0, (x,y,_)->(x,y), domain, false)
NoBody3D(domain) = Body((x,y,z,_=0)->-1.0, (x,y,z,_)->(x,y,z), domain, false)

Base.:+(a::Body, b::Body) = Body((x,t)->min(a.sdf(x,t), b.sdf(x,t)), (x,t)->ifelse(a.sdf(x,t) < b.sdf(x,t), a.map(x,t), b.map(x,t)), a.domain, a.compose || b.compose)
∩(a::Body, b::Body) = Body((x,t)->max(a.sdf(x,t), b.sdf(x,t)), (x,t)->ifelse(a.sdf(x,t) > b.sdf(x,t), a.map(x,t), b.map(x,t)), a.domain, a.compose || b.compose)
⊖(a::Body, b::Body) = Body((x,t)->max(a.sdf(x,t), -b.sdf(x,t)), a.map, a.domain, a.compose)
c(a::Body) = Body((x,t)->-a.sdf(x,t), a.map, a.domain, false)

measure(body::Body, x, t) = measure(body.sdf, body.map, x, t)
curvature(body::Body, point) = curvature(body.sdf, point)

"""
    d, n, v = measure(sdf, map, x, t)

Determine the signed distance, normal vector, and velocity of a body at a given point.
"""
function measure(sdf,map,x,t)

    # Determine the signed distance
    d = sdf(x...,t)

    # Determine the normal vector
    n = ForwardDiff.gradient(x -> sdf(x...,t), x)
    m = norm(n)
    n = n ./ m

    # Determine the velocity
    
    # Convert the map from a function of scalar x, y and t to a function of vector x and t
    #map = (x,t) -> map(x...,t)
    #sdf = (x,t) -> sdf(x...,t)

    # Compute the velocity
    J = ForwardDiff.jacobian(x -> map(x,t), x)
    dot = ForwardDiff.derivative(t -> map(x,t), t)
    v = -J\dot

    return (d,n,v)
end

function curvature(sdf, point)
    # Compute the gradient of the signed distance function at the point
    grad = ForwardDiff.gradient(x -> sdf(x...), point)

    # Compute the Hessian of the signed distance function at the point
    hess = ForwardDiff.hessian(x -> sdf(x...), point)
    adjoint_hess = det(hess) * inv(hess)

    # Compute the gaussian curvature
    gaussian_curvature = grad' * adjoint_hess * grad / norm(grad)^4

    # Compute the mean curvature
    mean_curvature = (grad' * hess * grad - norm(grad)^2 * tr(hess)) / (2*norm(grad)^3)

    return gaussian_curvature, mean_curvature
end
