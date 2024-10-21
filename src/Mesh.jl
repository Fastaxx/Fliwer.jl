abstract type AbstractMesh end

# Define the CartesianMesh structure
struct CartesianMesh{N}<: AbstractMesh
    h::NTuple{N, Array{Float64}}
    x0::NTuple{N, Float64}

    # Constructor for a random Cartesian mesh
    CartesianMesh(h::NTuple{N, Array{Float64}}, x0::NTuple{N, Float64}) where N = new{N}(h, x0)

    # Constructor for a uniform Cartesian mesh
    function CartesianMesh(N::Int, n::Array{Int}, L::Array{Float64}, x0::Array{Float64})
        h = [L[i] / n[i] for i in 1:N]
        return new{N}(h, x0)
    end
end
nC(mesh::CartesianMesh{N}) where N = prod(length.(mesh.h))


function centers(mesh::CartesianMesh{N}) where N
    center = [cumsum([mesh.x0[i]; mesh.h[i]]) for i in 1:N]
    return tuple(center...)
end

function nodes(mesh::CartesianMesh{N}) where N
    center = [cumsum([mesh.x0[i]; mesh.h[i]]) for i in 1:N]
    node = [(center[i][1:end-1] + center[i][2:end]) / 2 for i in 1:N]
    return tuple(node...)
end

"""
function edges(mesh::CartesianMesh{N}) where N
    nodes = [cumsum([mesh.x0[i]; mesh.h[i]]) for i in 1:N]
    edges = [(nodes[i][1:end-1] + nodes[i][2:end]) / 2 for i in 1:N]
    return tuple(edges...)
end

function faces(mesh::CartesianMesh{N}) where N
    nodes = [cumsum([mesh.x0[i]; mesh.h[i]]) for i in 1:N]
    faces = [(nodes[i][1:end-1] + nodes[i][2:end]) / 2 for i in 1:N]
    return tuple(faces...)
end
"""