abstract type AbstractMesh end

# Mesh struct 
# Type of mesh : x--!--x--!--x with x: cell center and !: cell boundary
# Start by a cell center then a cell boundary and so on  ... and finish by a cell center
struct Mesh{N} <: AbstractMesh
    nodes::NTuple{N, Vector{Float64}}
    centers::NTuple{N, Vector{Float64}}
    sizes::NTuple{N, Vector{Float64}}

    function Mesh(nodes::NTuple{N, AbstractVector{<:Real}}) where N
        centers = ntuple(i -> collect(Float64.(nodes[i])), N)
        nodes = ntuple(i -> (centers[i][1:end-1] .+ centers[i][2:end]) ./ 2.0, N)
        # Compute the sizes of the cells : The first and last cells have a size that is half of the others
        sizes = ntuple(i -> diff(nodes[i]), N)
        sizes = ntuple(i -> [sizes[i][1] / 2; sizes[i][1:end]; sizes[i][end] / 2], N)
        return new{N}(nodes, centers, sizes)
    end
end

# Function to extract border cells from a Mesh
function get_border_cells(mesh::Mesh{N}) where N
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
nC(mesh::Mesh{N}) where N = prod(length.(mesh.centers))

# --- Examples below ---

# For a 1D mesh:
x = collect(range(0.0, stop=1.0, length=5))
mesh1D = Mesh((x,))
borders1D = get_border_cells(mesh1D)
@show mesh1D.centers
@show nC(mesh1D)
@show mesh1D.sizes
@show length(mesh1D.sizes[1])

# For a 2D mesh:
x = collect(range(0.0, stop=1.0, length=11))
y = collect(range(0.0, stop=1.0, length=11))
mesh2D = Mesh((x, y))
borders2D = get_border_cells(mesh2D)
#@show borders2D
@show mesh2D.sizes

# For a 3D mesh:
x = collect(range(0.0, stop=1.0, length=11))
y = collect(range(0.0, stop=1.0, length=11))
z = collect(range(0.0, stop=1.0, length=11))
mesh3D = Mesh((x, y, z))
borders3D = get_border_cells(mesh3D)
#@show borders3D

# For a 4D mesh (if needed):
x = collect(range(0.0, stop=1.0, length=11))
y = collect(range(0.0, stop=1.0, length=11))
z = collect(range(0.0, stop=1.0, length=11))
t = collect(range(0.0, stop=1.0, length=11))
mesh4D = Mesh((x, y, z, t))
borders4D = get_border_cells(mesh4D)
#@show borders4D

# Define the body with a signed distance function :
# Two ways to define the same function (vectorized or not) and convert one to the other
# VOFI uses the non-vectorized version and ImplicitIntegration uses the vectorized version
Φ(X) = sqrt(X[1]^2 + X[2]^2) - 0.5
ϕ(x, y) = Φ([x, y])

LS(x,y) = sqrt(x^2 + y^2) - 0.5
ls(X) = LS(X[1], X[2])

@show ϕ(0.5, 0.5), Φ([0.5, 0.5])
@show ls([0.5, 0.5]), LS(0.5, 0.5)

# Define a moving body : a circle moving along the x-axis : vectorized version
ϕ(X) = sqrt((X[1] - 0.5)^2 + (X[2] - 0.5)^2) - 0.2
map(X,t) = [X[1] + t, X[2]]
Φ(X,t) = ϕ(map(X,t))

@show Φ([0.5, 0.5], 0.0), Φ([0.5, 0.5], 0.1)

# Define a moving body : a circle moving along the x-axis : non-vectorized version
ϕ(x,y) = sqrt((x - 0.5)^2 + (y - 0.5)^2) - 0.2
map(x,y,t) = x + t, y
Φ(x,y,t) = ϕ(map(x,y,t))

@show Φ(0.5, 0.5, 0.0), Φ(0.5, 0.5, 0.1)