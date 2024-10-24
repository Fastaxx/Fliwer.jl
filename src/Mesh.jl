abstract type AbstractMesh end


"""
struct Mesh{N}

The `Mesh` struct represents a mesh in N-dimensional space.
x--!--x--!--x--!--x--!--x--!--x
x: cell center
!: cell boundary

# Fields
- `h::NTuple{N, Vector{Float64}}`: Sizes of the cells for each dimension.
- `x0::NTuple{N, Float64}`: Origin for each dimension.
- `nodes::NTuple{N, Vector{Float64}}`: Nodes (cell boundaries) for each dimension.
- `centers::NTuple{N, Vector{Float64}}`: Centers of the cells for each dimension.
"""
struct CartesianMesh{N} <: AbstractMesh
    h::NTuple{N, Vector{Float64}}        # Cell sizes for each dimension
    x0::NTuple{N, Float64}               # Origin for each dimension
    nodes::NTuple{N, Vector{Float64}}    # Nodes (cell edges) for each dimension
    centers::NTuple{N, Vector{Float64}}  # Cell centres for each dimension

    # Principal constructor
    function CartesianMesh(h::NTuple{N, Vector{Float64}}, x0::NTuple{N, Float64}) where N
        center = ntuple(i -> cumsum([x0[i]; h[i]]), N)
        nodes = ntuple(i -> (center[i][1:end-1] .+ center[i][2:end]) ./ 2.0, N)

        return new{N}(h, x0, nodes, center)
    end

     # Constructeur pour un maillage cartésien uniforme avec moins de paramètres
    function CartesianMesh(n::NTuple{N, Int}, domain_size::NTuple{N, Float64}, 
        x0::NTuple{N, Float64}=ntuple(_ -> 0.0, N)) where N
        # Calcul des tailles des cellules uniformes
        h_uniform = ntuple(i -> fill(domain_size[i] / n[i], n[i]), N)

        # Calcul des centres des cellules
        centers_uniform = ntuple(i -> [x0[i] + (j + 0.5) * (domain_size[i] / n[i]) for j in 0:(n[i]-1)], N)

        # Calcul des noeuds 
        nodes_uniform  = ntuple(i -> [x0[i] + j * (domain_size[i] / n[i]) for j in 0:n[i]], N)

        # Création de l'instance en utilisant l'constructeur interne
        new{N}(h_uniform, x0, nodes_uniform, centers_uniform)
    end
end

nC(mesh::CartesianMesh{N}) where N = prod(length.(mesh.h))

