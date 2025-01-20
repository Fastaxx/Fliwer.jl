"""
    struct MeshTag

A mutable struct representing a mesh tag.

# Fields
- `border_cells`: An array of tuples containing the Cartesian indices and corresponding integers for border cells.
- `cut_cells`: An array of tuples containing the Cartesian indices and corresponding integers for cut cells.
- `regular_cells`: An array of tuples containing the Cartesian indices and corresponding integers for regular cells.
"""
mutable struct MeshTag
    border_cells::Array{Tuple{CartesianIndex, Int}, 1}
    cut_cells::Array{Tuple{CartesianIndex, Int}, 1}
    regular_cells::Array{Tuple{CartesianIndex, Int}, 1}
end

"""
    abstract type AbstractMesh

An abstract type representing a mesh.
"""
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
- `faces::NTuple{N, NTuple{N, Vector{Float64}}}`: Faces for each dimension.
- `tag::MeshTag`: Tagging information for cells.
"""
mutable struct CartesianMesh{N} <: AbstractMesh
    h::NTuple{N, Vector{Float64}}        # Cell sizes for each dimension
    x0::NTuple{N, Float64}               # Origin for each dimension
    nodes::NTuple{N, Vector{Float64}}    # Nodes (cell edges) for each dimension
    centers::NTuple{N, Vector{Float64}}  # Cell centres for each dimension
    faces::NTuple{N, NTuple{N, Vector{Float64}}}    # Faces for each dimension
    tag::MeshTag                         # Tagging information for cells

    # Main constructor
    function CartesianMesh(h::NTuple{N, Vector{Float64}}, x0::NTuple{N, Float64}) where N
        println("=== Creating CartesianMesh ===\n- h: $(h)\n- x0: $(x0)\n")

        center = ntuple(i -> cumsum([x0[i]; h[i]]), N)
        nodes = ntuple(i -> (center[i][1:end-1] .+ center[i][2:end]) ./ 2.0, N)
        
        # Calcul des faces
        # faces = (facex, facey, facez, ...)
        # Chaque facex, facey, etc., est un tuple de positions pour chaque dimension
        faces = ntuple(dim -> ntuple(i -> (dim == i) ? 
            [v + 0.5 * h[i][j] for (j, v) in enumerate(center[i][1:end-1])]
            : center[i], N), N)

        tag = MeshTag([], [], [])

        return new{N}(h, x0, nodes, center, faces, tag)
    end

     # Constructor for a uniform Cartesian mesh with fewer parameters
    function CartesianMesh(n::NTuple{N, Int}, domain_size::NTuple{N, Float64}, 
        x0::NTuple{N, Float64}=ntuple(_ -> 0.0, N)) where N
        println("=== Creating uniform CartesianMesh ===\n- n: $(n)\n- domain_size: $(domain_size)\n- x0: $(x0)\n")

        h_uniform = ntuple(i -> fill(domain_size[i] / n[i], n[i]), N)
        centers_uniform = ntuple(i -> [x0[i] + j * (domain_size[i] / n[i]) for j in 0:n[i]-1], N)
        nodes_uniform  = ntuple(i -> [x0[i] + (j + 0.5) * (domain_size[i] / n[i]) for j in 0:(n[i])], N) 
        # Calcul des faces
        faces_uniform = ntuple(dim -> ntuple(i -> (dim == i) ? 
            [centers_uniform[i][j] + 0.5 * (domain_size[i] / n[i]) for j in 1:n[i]]
            : centers_uniform[i], N), N)

        tag = MeshTag([], [], [])

        new{N}(h_uniform, x0, nodes_uniform, centers_uniform, faces_uniform, tag)
    end
end

nC(mesh::CartesianMesh{N}) where N = prod(length.(mesh.h))

"""
    mutable struct CartesianSpaceTimeMesh{N} <: AbstractMesh

A simple (N+1)-dimensional mesh that combines an existing CartesianMesh{N}
(covering only space) with a discrete time vector.
Fields mirror the structure of a usual `CartesianMesh` but add time as the (N+1)-th dimension.

Use the provided constructor `CartesianSpaceTimeMesh(spaceMesh::CartesianMesh{N}, t::Vector{Float64})`
to build one from a spatial mesh and a list of time instants.
"""
mutable struct CartesianSpaceTimeMesh{M} <: AbstractMesh
    h         :: NTuple{M, Vector{Float64}}
    x0        :: NTuple{M, Float64}
    nodes     :: NTuple{M, Vector{Float64}}
    centers   :: NTuple{M, Vector{Float64}}
    faces     :: NTuple{M, NTuple{M, Vector{Float64}}}
    tag       :: MeshTag

    function CartesianSpaceTimeMesh(spaceMesh::CartesianMesh{N}, t::Vector{Float64}; tag = MeshTag([], [], [])) where {N}
        # Compute M = N + 1 at runtime
        local M = N + 1

        # -- Build time arrays --
        Δt          = [t[i+1] - t[i] for i in 1:length(t)-1]
        x0_time     = t[1]
        centers_time = [(t[i+1] + t[i]) / 2 for i in 1:length(t)-1]

        # -- Assemble the ND+1 arrays --
        new_h       = ntuple(i -> i <= N ? spaceMesh.h[i] : Δt, M)
        new_x0      = ntuple(i -> i <= N ? spaceMesh.x0[i] : x0_time, M)
        new_nodes   = ntuple(i -> i <= N ? spaceMesh.nodes[i] : t, M)
        new_centers = ntuple(i -> i <= N ? spaceMesh.centers[i] : centers_time, M)
        new_faces   = ntuple(dim -> ntuple(i -> (dim == i) ? 
            [new_centers[i][j] + 0.5 * new_h[i][j] for j in 1:length(new_h[i])]
            : new_centers[i], M), M)

        return new{M}(new_h, new_x0, new_nodes, new_centers, new_faces, tag)
    end
end

nC(stmesh::CartesianSpaceTimeMesh{N}) where N = prod(length.(stmesh.h))
