abstract type AbstractCapacity end

mutable struct Capacity{N} <: AbstractCapacity
    A :: NTuple{N, SparseMatrixCSC{Float64, Int}}
    B :: NTuple{N, SparseMatrixCSC{Float64, Int}}
    V :: SparseMatrixCSC{Float64, Int}
    W :: NTuple{N, SparseMatrixCSC{Float64, Int}}
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
    return A, B, V, W
end

function Capacity(body::AbstractBody, mesh::CartesianMesh)
    A, B, V, W = VOFI(body, mesh)
    
    N = length(A)
    
    return Capacity{N}(A, B, V, W, mesh)
end

function measure!(capacity::AbstractCapacity, body::AbstractBody; t=0)
    A, B, V, W = VOFI(body, capacity.mesh)
    capacity.A, capacity.B, capacity.V, capacity.W = A, B, V, W
end