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

Base.:+(a::Body, b::Body) = Body((x,t)->min(a.sdf(x,t), b.sdf(x,t)), (x,t)->ifelse(a.sdf(x,t) < b.sdf(x,t), a.map(x,t), b.map(x,t)), a.domain, a.compose || b.compose)
∩(a::Body, b::Body) = Body((x,t)->max(a.sdf(x,t), b.sdf(x,t)), (x,t)->ifelse(a.sdf(x,t) > b.sdf(x,t), a.map(x,t), b.map(x,t)), a.domain, a.compose || b.compose)
⊖(a::Body, b::Body) = Body((x,t)->max(a.sdf(x,t), -b.sdf(x,t)), a.map, a.domain, a.compose)
c(a::Body) = Body((x,t)->-a.sdf(x,t), a.map, a.domain, false)

function measure(sdf, map, x, t)
    # Définir une fonction qui ne dépend que de x, en fixant t
    f(x_vec) = sdf(x_vec, t)
    
    # Calculer la distance signée
    d = f(x)
    
    # Calculer le gradient de la distance signée par rapport à x
    ∇d = ForwardDiff.gradient(f, x)
    
    # Vérifier si le gradient est nul pour éviter la division par zéro
    if norm(∇d) == 0
        n = zeros(Float64, length(x))
    else
        # Calculer le vecteur normal
        n = ∇d ./ norm(∇d)
    end
    
    # Calculer la matrice jacobienne de la fonction map par rapport à x
    J = ForwardDiff.jacobian(x_vec -> map(x_vec, t), x)
    
    # Calculer la dérivée temporelle de la fonction map
    dot = ForwardDiff.derivative(t_val -> map(x, t_val), t)
    
    # Calculer la vitesse
    v = -J \ dot
    
    return d, n, v
end
