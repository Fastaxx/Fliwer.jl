using SparseArrays

abstract type AbstractBody end

struct NoBody <: AbstractBody end
abstract type AbstractFluid end

# Routine pour évaluer les capacités
function evaluate_capacity(body::AbstractBody, ::Type{T}, ::Val{N}) where {N, T}
    # Exemple de matrices pour 2D
    A = ntuple(i -> sparse([1.0 0.0; 0.0 1.0]), N)
    B = ntuple(i -> sparse([1.0 0.0; 0.0 1.0]), N)
    V = sparse([1.0 0.0; 0.0 1.0])
    W = ntuple(i -> sparse([1.0 0.0; 0.0 1.0]), N)
    return A, B, V, W
end

struct Fluid{N,T} <: AbstractFluid
    body::AbstractBody
    A::NTuple{N,SparseMatrixCSC{T,Int}}
    B::NTuple{N,SparseMatrixCSC{T,Int}}
    V::SparseMatrixCSC{T,Int}
    W::NTuple{N,SparseMatrixCSC{T,Int}}
end

# Constructeur externe pour Fluid
function Fluid{N,T}(body::AbstractBody) where {N,T}
    A, B, V, W = evaluate_capacity(body, T, Val(N))
    Fluid{N,T}(body, A, B, V, W)
end

# Test
body = NoBody()
fluid = Fluid{2,Float64}(body)

println(fluid)

println(fluid.body)
println(fluid.A)
println(fluid.B)
println(fluid.V)
println(fluid.W)