abstract type AbstractCapacity end

mutable struct Capacity{N} <: AbstractCapacity
    A :: NTuple{N, SparseMatrixCSC{Float64, Int}}   # A capacity : Ax, Ay
    B :: NTuple{N, SparseMatrixCSC{Float64, Int}}   # B capacity : Bx, By
    V :: SparseMatrixCSC{Float64, Int}              # Volume
    W :: NTuple{N, SparseMatrixCSC{Float64, Int}}   # Staggered Volume
    C :: Vector{SVector{N,Float64}}                 # Cell Centroid
    mesh :: CartesianMesh{N}
end

function VOFI(body::AbstractBody, mesh::CartesianMesh)
    N = length(mesh.h)
    nc = nC(mesh)

    Vs, bary = spzeros(nc), zeros(N)
    As, Bs, Ws = (spzeros(nc), spzeros(nc), spzeros(nc)), (spzeros(nc), spzeros(nc), spzeros(nc)), (spzeros(nc), spzeros(nc), spzeros(nc))

    Vs, bary = integrate(Tuple{0}, body.sdf, mesh.centers, Float64, zero)
    As = integrate(Tuple{1}, body.sdf, mesh.centers, Float64, zero)
    Ws = integrate(Tuple{0}, body.sdf, mesh.centers, Float64, zero, bary)
    Bs = integrate(Tuple{1}, body.sdf, mesh.centers, Float64, zero, bary)

    C = bary
    if N == 1
        V = spdiagm(0 => Vs)
        A = (spdiagm(0 => As[1]),)
        B = (spdiagm(0 => Bs[1]),)
        W = (spdiagm(0 => Ws[1]),)
    elseif N == 2
        V = spdiagm(0 => Vs)
        A = (spdiagm(0 => As[1]), spdiagm(0 => As[2]))
        B = (spdiagm(0 => Bs[1]), spdiagm(0 => Bs[2]))
        W = (spdiagm(0 => Ws[1]), spdiagm(0 => Ws[2]))
    elseif N == 3
        V = spdiagm(0 => Vs)
        A = (spdiagm(0 => As[1]), spdiagm(0 => As[2]), spdiagm(0 => As[3]))
        B = (spdiagm(0 => Bs[1]), spdiagm(0 => Bs[2]), spdiagm(0 => Bs[3]))
        W = (spdiagm(0 => Ws[1]), spdiagm(0 => Ws[2]), spdiagm(0 => Ws[3]))
    end
    return A, B, V, W, C
end

function Capacity(body::AbstractBody, mesh::CartesianMesh)
    A, B, V, W, C = VOFI(body, mesh)
    
    N = length(A)
    
    return Capacity{N}(A, B, V, W, C, mesh)
end

function measure!(capacity::AbstractCapacity, body::AbstractBody; t=0)
    A, B, V, W, C = VOFI(body, capacity.mesh)
    capacity.A, capacity.B, capacity.V, capacity.W, capacity.C = A, B, V, W, C
end