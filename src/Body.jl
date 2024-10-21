"""
    AbstractBody

Embedded body Abstract type. Any `AbstractBody` subtype must implement

d = sdf(body::AbstractBody, x, t=0)

and

d,n,V = measure(body::AbstractBody, x, t=0)

where `d` is the signed distance from `x` to the body at time `t`,
and `n` & `V` are the normal and velocity vectors implied at `x`.
"""
abstract type AbstractBody end

"""
    NoBody

Use for a simulation without any body.
"""
struct NoBody <: AbstractBody end

"""
    measure!(fluid::AbstractFluid, body::AbstractBody; t=0)

Queries the body geometry to fill the arrays:

- `fluid.A`, A geometric capacity
- `fluid.B`, B geometric capacity
- `fluid.V`, V geometric capacity
- `fluid.W`, W geometric capacity

at time `t`.
"""
