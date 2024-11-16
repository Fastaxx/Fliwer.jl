"""
    abstract type AbstractCapacity

Abstract type representing a capacity.
"""
abstract type AbstractCapacity end

"""
    mutable struct Capacity{N} <: AbstractCapacity

The `Capacity` struct represents the capacity of a system in `N` dimensions.

# Fields
- `A`: A capacity represented by `N` sparse matrices (`Ax`, `Ay`).
- `B`: B capacity represented by `N` sparse matrices (`Bx`, `By`).
- `V`: Volume capacity represented by a sparse matrix.
- `W`: Staggered volume capacity represented by `N` sparse matrices.
- `C_ω`: Cell centroid represented by a vector of `N`-dimensional static vectors.
- `C_γ`: Interface centroid represented by a vector of `N`-dimensional static vectors.
- `Γ`: Interface norm represented by a sparse matrix.
- `mesh`: Cartesian mesh of `N` dimensions.

"""
mutable struct Capacity{N} <: AbstractCapacity
    A :: NTuple{N, SparseMatrixCSC{Float64, Int}}   # A capacity : Ax, Ay
    B :: NTuple{N, SparseMatrixCSC{Float64, Int}}   # B capacity : Bx, By
    V :: SparseMatrixCSC{Float64, Int}              # Volume
    W :: NTuple{N, SparseMatrixCSC{Float64, Int}}   # Staggered Volume
    C_ω :: Vector{SVector{N,Float64}}               # Cell Centroid
    C_γ :: Vector{SVector{N,Float64}}               # Interface Centroid
    Γ :: SparseMatrixCSC{Float64, Int}              # Interface Norm
    mesh :: CartesianMesh{N}
end

"""
    VOFI(body::AbstractBody, mesh::CartesianMesh)

Compute the Capacity quantities based on VOFI for a given body and mesh.

# Arguments
- `body::AbstractBody`: The body for which to compute the VOFI quantities.
- `mesh::CartesianMesh`: The mesh on which to compute the VOFI quantities.

# Returns
- `A::Tuple`: The A matrices for each dimension of the mesh.
- `B::Tuple`: The B matrices for each dimension of the mesh.
- `V::SparseMatrixCSC`: The V matrix.
- `W::Tuple`: The W matrices for each dimension of the mesh.
- `C_ω::Vector`: The C_ω vector : Cell centroid.
- `C_γ::Vector`: The C_γ vector : Interface centroid.
- `Γ::SparseMatrixCSC`: The Γ matrix : Interface Norm.

"""
function VOFI(body::AbstractBody, mesh::CartesianMesh)
    N = length(mesh.h)
    nc = nC(mesh)

    Vs, bary, interface_length = spzeros(nc), zeros(N), spzeros(nc)
    As, Bs, Ws = (spzeros(nc), spzeros(nc), spzeros(nc)), (spzeros(nc), spzeros(nc), spzeros(nc)), (spzeros(nc), spzeros(nc), spzeros(nc))

    Vs, bary = integrate(Tuple{0}, body.sdf, mesh.nodes, Float64, zero)
    As = integrate(Tuple{1}, body.sdf, mesh.nodes, Float64, zero)
    Ws = integrate(Tuple{0}, body.sdf, mesh.nodes, Float64, zero, bary)
    Bs = integrate(Tuple{1}, body.sdf, mesh.nodes, Float64, zero, bary)

    C_ω = bary
    if N == 1
        V = spdiagm(0 => Vs)
        A = (spdiagm(0 => As[1]),)
        B = (spdiagm(0 => Bs[1]),)
        W = (spdiagm(0 => Ws[1]),)
        Γ = spdiagm(0 => interface_length)
    elseif N == 2
        V = spdiagm(0 => Vs)
        A = (spdiagm(0 => As[1]), spdiagm(0 => As[2]))
        B = (spdiagm(0 => Bs[1]), spdiagm(0 => Bs[2]))
        W = (spdiagm(0 => Ws[1]), spdiagm(0 => Ws[2]))
        Γ = spdiagm(0 => interface_length)
    elseif N == 3
        V = spdiagm(0 => Vs)
        A = (spdiagm(0 => As[1]), spdiagm(0 => As[2]), spdiagm(0 => As[3]))
        B = (spdiagm(0 => Bs[1]), spdiagm(0 => Bs[2]), spdiagm(0 => Bs[3]))
        W = (spdiagm(0 => Ws[1]), spdiagm(0 => Ws[2]), spdiagm(0 => Ws[3]))
        Γ = spdiagm(0 => interface_length)
    end

    C_γ = compute_c_γ(mesh.nodes, body.sdf)
    return A, B, V, W, C_ω, C_γ, Γ
end

"""
    Capacity(body::AbstractBody, mesh::CartesianMesh)

Compute the capacity of a body in a given mesh using the VOFI method.

# Arguments
- `body::AbstractBody`: The body for which to compute the capacity.
- `mesh::CartesianMesh`: The mesh in which the body is located.

# Returns
- `Capacity{N}`: The capacity of the body.
"""
function Capacity(body::AbstractBody, mesh::CartesianMesh)
    A, B, V, W, C_ω, C_γ, Γ = VOFI(body, mesh)
    
    N = length(A)
    
    return Capacity{N}(A, B, V, W, C_ω, C_γ, Γ, mesh)
end

function measure!(capacity::AbstractCapacity, body::AbstractBody; t=0)
    A, B, V, W, C_ω, C_γ, Γ = VOFI(body, capacity.mesh)
    capacity.A, capacity.B, capacity.V, capacity.W, capacity.C_ω, capacity.C_γ, capacity.Γ = A, B, V, W, C_ω, C_γ, Γ
end