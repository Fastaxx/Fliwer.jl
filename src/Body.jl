
"""
    AbstractBody

Embedded body Abstract type. Any `AbstractBody` subtype must implement

d = sdf(body::AbstractBody, x, t=0)

and

d,n,V = measure(body::AbstractBody, x, t=0, fastd²=Inf)

where `d` is the signed distance from `x` to the body at time `t`,
and `n` & `V` are the normal and velocity vectors implied at `x`.
A fast-approximate method can return `≈d,zero(x),zero(x)` if `d^2>fastd²`.
"""
abstract type AbstractBody end

"""
    NoBody

Use for a simulation without any body.
"""
struct NoBody <: AbstractBody end

"""
    measure_sdf!(a::AbstractArray, body::AbstractBody, t=0)

Fill the array `a` with the signed distance field from `body` at time `t`.
"""
measure_sdf!(a::AbstractArray,body::AbstractBody,t=0;kwargs...) = @inside a[I] = sdf(body,loc(0,I,eltype(a)),t;kwargs...)

