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

function create_staggered_meshes(mesh_cc::CartesianMesh{1})
    # Extract cell centers and sizes
    x_centers = mesh_cc.centers[1]
    h_x = mesh_cc.h[1]

    # Compute staggered nodes for x-velocity (mesh_stx)
    xx = vcat(
        x_centers[1] - h_x[1]/2,
        x_centers[1:end-1] .+ h_x[1:end-1]/2,
        x_centers[end] + h_x[end]/2
    )

    # Create mesh_stx
    nx_stx = length(xx) - 1
    mesh_stx = CartesianMesh(
        (nx_stx,),
        (xx[end] - xx[1],),
        (xx[1],)
    )

    return (mesh_cc, mesh_stx)
end

function create_staggered_meshes(mesh_cc::CartesianMesh{2})
    # Extract cell centers and sizes
    x_centers = mesh_cc.centers[1]
    y_centers = mesh_cc.centers[2]
    h_x = mesh_cc.h[1]
    h_y = mesh_cc.h[2]

    # Compute staggered nodes for x-velocity (mesh_stx)
    xx = vcat(
        x_centers[1] - h_x[1]/2,
        x_centers[1:end-1] .+ h_x[1:end-1]/2,
        x_centers[end] + h_x[end]/2
    )

    # Create mesh_stx
    nx_stx = length(xx) - 1
    ny_stx = length(y_centers)
    mesh_stx = CartesianMesh(
        (nx_stx, ny_stx),
        (xx[end] - xx[1], y_centers[end] - y_centers[1]),
        (xx[1], y_centers[1])
    )

    # Compute staggered nodes for y-velocity (mesh_sty)
    yy = vcat(
        y_centers[1] - h_y[1]/2,
        y_centers[1:end-1] .+ h_y[1:end-1]/2,
        y_centers[end] + h_y[end]/2
    )

    # Create mesh_sty
    nx_sty = length(x_centers)
    ny_sty = length(yy) - 1
    mesh_sty = CartesianMesh(
        (nx_sty, ny_sty),
        (x_centers[end] - x_centers[1], yy[end] - yy[1]),
        (x_centers[1], yy[1])
    )

    return (mesh_cc, mesh_stx, mesh_sty)
end