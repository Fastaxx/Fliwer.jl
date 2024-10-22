abstract type AbstractFluid end

# Définition de la structure Fluid paramétrée par N et fixant T à Float64
mutable struct Fluid{N} <: AbstractFluid
    A :: NTuple{N, SparseMatrixCSC{Float64, Int}}
    B :: NTuple{N, SparseMatrixCSC{Float64, Int}}
    V :: SparseMatrixCSC{Float64, Int}
    W :: NTuple{N, SparseMatrixCSC{Float64, Int}}
    mesh :: CartesianMesh{N}
end

# Définition de la fonction VOFI qui retourne des NTuple
function VOFI(body::AbstractBody, mesh::CartesianMesh)
    # Déduction de la dimension N à partir du maillage
    N = length(mesh.h)  # Puisque h est un NTuple{N, Array{Float64}}
    nc = nC(mesh)
    mesh = centers(mesh)

    Vs, bary = spzeros(nc), zeros(N)
    As, Bs, Ws = (spzeros(nc), spzeros(nc), spzeros(nc)), (spzeros(nc), spzeros(nc), spzeros(nc)), (spzeros(nc), spzeros(nc), spzeros(nc))

    Vs, bary = integrate(Tuple{0}, body.sdf, mesh, Float64, zero)
    As = integrate(Tuple{1}, body.sdf, mesh, Float64, zero)
    Ws = integrate(Tuple{0}, body.sdf, mesh, Float64, zero, bary)
    Bs = integrate(Tuple{1}, body.sdf, mesh, Float64, zero, bary)
    
    # Génération des tuples A, B, W avec N matrices creuses
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

# Définition du constructeur extérieur pour Fluid
function Fluid(body::AbstractBody, mesh::CartesianMesh)
    A, B, V, W = VOFI(body, mesh)
    
    # Déduction de la dimension N à partir de la longueur du tuple A
    N = length(A)  # Puisque A est un NTuple{N, ...}
    
    # Création de l'instance de Fluid{N}
    return Fluid{N}(A, B, V, W, mesh)
end

"""
    measure!(fluid::AbstractFluid, body::AbstractBody; t=0)

Queries the body geometry to fill the arrays:

- `fluid.A`, A geometric capacity
- `fluid.B`, B geometric capacity
- `fluid.V`, V geometric capacity
- `fluid.W`, W geometric capacity

at time `t`.
"""
function measure!(fluid::AbstractFluid, body::AbstractBody; t=0)
    # Appeler VOFI pour obtenir les nouvelles capacités
    A_new, B_new, V_new, W_new = VOFI(body, fluid.mesh)
    
    # Mettre à jour les capacités de fluid
    fluid.A = A_new
    fluid.B = B_new
    fluid.V = V_new
    fluid.W = W_new
    
    return nothing
end