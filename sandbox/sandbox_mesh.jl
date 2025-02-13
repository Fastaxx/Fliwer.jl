# ScalarMesh struct 
# Type of mesh : x--!--x--!--x with x: cell center and !: cell boundary
# Start by a cell center then a cell boundary and so on  ... and finish by a cell center
struct ScalarMesh{N}
    nodes::NTuple{N, Vector{Float64}}
    centers::NTuple{N, Vector{Float64}}
    sizes::NTuple{N, Vector{Float64}}

    function ScalarMesh(nodes::NTuple{N, AbstractVector{<:Real}}) where N
        centers = ntuple(i -> collect(Float64.(nodes[i])), N)
        nodes = ntuple(i -> (centers[i][1:end-1] .+ centers[i][2:end]) ./ 2.0, N)
        # Compute the sizes of the cells : The first and last cells have a size that is half of the others
        sizes = ntuple(i -> diff(nodes[i]), N)
        sizes = ntuple(i -> [sizes[i][1] / 2; sizes[i][1:end]; sizes[i][end] / 2], N)
        return new{N}(nodes, centers, sizes)
    end
end

# Function to extract border cells from a ScalarMesh
function get_border_cells(mesh::ScalarMesh{N}) where N
    # Number of cells in each dimension
    dims = ntuple(i -> length(mesh.centers[i]), N)
    border_cells = Vector{Tuple{CartesianIndex, NTuple{N, Float64}}}()
    
    # Iterate over all cell indices using Iterators.product
    for idx in Iterators.product((1:d for d in dims)...)
        # A cell is at the border if any index equals 1 or the maximum in that dimension.
        if any(d -> idx[d] == 1 || idx[d] == dims[d], 1:N)
            # Get the physical cell center: tuple (mesh.centers[1][i₁], mesh.centers[2][i₂], ...)
            pos = ntuple(d -> mesh.centers[d][idx[d]], N)
            push!(border_cells, (CartesianIndex(idx), pos))
        end
    end
    return border_cells
end

# Function to get the total number of cells in a mesh
nC(mesh::ScalarMesh{N}) where N = prod(length.(mesh.centers))

# --- Examples below ---

# For a 1D mesh:
x = collect(range(0.0, stop=1.0, length=5))
mesh1D = ScalarMesh((x,))
borders1D = get_border_cells(mesh1D)
@show mesh1D.centers
@show nC(mesh1D)
@show mesh1D.sizes
@show length(mesh1D.sizes[1])

# For a 2D mesh:
x = collect(range(0.0, stop=1.0, length=11))
y = collect(range(0.0, stop=1.0, length=11))
mesh2D = ScalarMesh((x, y))
borders2D = get_border_cells(mesh2D)
#@show borders2D
@show mesh2D.sizes

# For a 3D mesh:
x = collect(range(0.0, stop=1.0, length=11))
y = collect(range(0.0, stop=1.0, length=11))
z = collect(range(0.0, stop=1.0, length=11))
mesh3D = ScalarMesh((x, y, z))
borders3D = get_border_cells(mesh3D)
#@show borders3D

# For a 4D mesh (if needed):
x = collect(range(0.0, stop=1.0, length=11))
y = collect(range(0.0, stop=1.0, length=11))
z = collect(range(0.0, stop=1.0, length=11))
t = collect(range(0.0, stop=1.0, length=11))
mesh4D = ScalarMesh((x, y, z, t))
borders4D = get_border_cells(mesh4D)
#@show borders4D
