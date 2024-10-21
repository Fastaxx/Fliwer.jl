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
    
    # Génération des tuples A, B, W avec N matrices creuses
    A = ntuple(i -> sparse(rand(nc, nc)), N)
    B = ntuple(i -> sparse(rand(nc, nc)), N)
    W = ntuple(i -> sparse(rand(nc, nc)), N)
    
    # Génération de la matrice V
    V = sparse(rand(nc, nc))
    
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
    
    # Mettre à jour les champs de Fluid
    fluid.A = A_new
    fluid.B = B_new
    fluid.V = V_new
    fluid.W = W_new
    
    return nothing
end