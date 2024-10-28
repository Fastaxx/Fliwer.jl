mutable struct MeshTag
    border_cells::Array{Tuple{CartesianIndex, Int}, 1}
    cut_cells::Array{Tuple{CartesianIndex, Int}, 1}
    regular_cells::Array{Tuple{CartesianIndex, Int}, 1}
end

function find_border(mesh)
    centers = mesh.centers
    dims = length(centers)
    border_cells = []

    if dims == 1
        for i in 1:length(centers[1])
            if i == 1 || i == length(centers[1])
                cartesian_index = CartesianIndex(i)
                linear_index = LinearIndices((length(centers[1]),))[cartesian_index]
                push!(border_cells, (cartesian_index, linear_index))
            end
        end
    elseif dims == 2
        for j in 1:length(centers[2])
            for i in 1:length(centers[1])
                if i == 1 || i == length(centers[1]) || j == 1 || j == length(centers[2])
                    cartesian_index = CartesianIndex(i, j)
                    linear_index = LinearIndices((length(centers[1]), length(centers[2])))[cartesian_index]
                    push!(border_cells, (cartesian_index, linear_index))
                end
            end
        end
    elseif dims == 3
        for k in 1:length(centers[3])
            for j in 1:length(centers[2])
                for i in 1:length(centers[1])
                    if i == 1 || i == length(centers[1]) || j == 1 || j == length(centers[2]) || k == 1 || k == length(centers[3])
                        cartesian_index = CartesianIndex(i, j, k)
                        linear_index = LinearIndices((length(centers[1]), length(centers[2]), length(centers[3])))[cartesian_index]
                        push!(border_cells, (cartesian_index, linear_index))
                    end
                end
            end
        end
    else
        error("Unsupported dimension: $dims")
    end

    return border_cells
end

function eval_sdf(mesh::CartesianMesh, body::Body)
    # Evaluate the signed distance function on the mesh.nodes
    nodes = mesh.nodes
    dim = length(nodes)
    sdf = zeros(prod(length.(nodes)))
    
    if dim == 1
        x = nodes[1]
        return [body.sdf(i) for i in x]
    elseif dim == 2
        x, y = nodes
        return [body.sdf(i, j) for i in x, j in y]
    elseif dim == 3
        x, y, z = nodes
        return [body.sdf(i, j, k) for i in x, j in y, k in z]
    else
        error("Unsupported dimension: $dim")
    end

end


# Use VOFI to identify the cut cells instead
function find_cut(mesh::CartesianMesh, body::Body)
    # Evaluate the signed distance function on the mesh
    sdf = eval_sdf(mesh, body)
    cut_cells = []
    dims = length(mesh.centers)
    # If the signed distance function changes sign in a cell, it is a cut cell
    if dims == 1
        for i in axes(sdf, 1)[begin:end-1]
            if sdf[i] * sdf[i+1] < 0
                cartesian_index = CartesianIndex(i)
                linear_index = LinearIndices((length(mesh.centers[1]),))[cartesian_index]
                push!(cut_cells, (cartesian_index, linear_index))
            end
        end
    elseif dims == 2
        for j in axes(sdf, 2)[begin:end-1]
            for i in axes(sdf, 1)[begin:end-1]
                if sdf[i, j] * sdf[i+1, j] < 0 ||
                    sdf[i, j] * sdf[i, j+1] < 0 ||
                    sdf[i+1, j] * sdf[i+1, j+1] < 0 ||
                    sdf[i, j+1] * sdf[i+1, j+1] < 0
                    cartesian_index = CartesianIndex(i, j)
                    linear_index = LinearIndices((length(mesh.centers[1]), length(mesh.centers[2])))[cartesian_index]
                    push!(cut_cells, (cartesian_index, linear_index))
                end
            end
        end
    elseif dims == 3
        for k in axes(sdf, 3)[begin:end-1]
            for j in axes(sdf, 2)[begin:end-1]
                for i in axes(sdf, 1)[begin:end-1]
                    if sdf[i, j, k] * sdf[i+1, j, k] < 0 ||
                        sdf[i, j, k] * sdf[i, j+1, k] < 0 ||
                        sdf[i, j, k] * sdf[i, j, k+1] < 0 ||
                        sdf[i+1, j, k] * sdf[i+1, j+1, k] < 0 ||
                        sdf[i+1, j, k] * sdf[i+1, j, k+1] < 0 ||
                        sdf[i, j+1, k] * sdf[i, j+1, k+1] < 0 ||
                        sdf[i, j, k+1] * sdf[i+1, j, k+1] < 0 ||
                        sdf[i, j, k+1] * sdf[i, j+1, k+1] < 0 ||
                        sdf[i+1, j+1, k] * sdf[i+1, j+1, k+1] < 0
                        cartesian_index = CartesianIndex(i, j, k)
                        linear_index = LinearIndices((length(mesh.centers[1]), length(mesh.centers[2]), length(mesh.centers[3])))[cartesian_index]
                        push!(cut_cells, (cartesian_index, linear_index))
                    end
                end
            end
        end
    else
        error("Unsupported dimension: $dims")
    end

    return cut_cells
end

# Use VOFI to identify the regular cells instead
function find_regular(mesh::CartesianMesh)
    centers = mesh.centers
    dims = length(centers)
    regular_cells = []

    if dims == 1
        for i in 2:length(centers[1])-1
            cartesian_index = CartesianIndex(i)
            linear_index = LinearIndices((length(centers[1]),))[cartesian_index]
            push!(regular_cells, (cartesian_index, linear_index))
        end
    elseif dims == 2
        for j in 2:length(centers[2])-1
            for i in 2:length(centers[1])-1
                cartesian_index = CartesianIndex(i, j)
                linear_index = LinearIndices((length(centers[1]), length(centers[2])))[cartesian_index]
                push!(regular_cells, (cartesian_index, linear_index))
            end
        end
    elseif dims == 3
        for k in 2:length(centers[3])-1
            for j in 2:length(centers[2])-1
                for i in 2:length(centers[1])-1
                    cartesian_index = CartesianIndex(i, j, k)
                    linear_index = LinearIndices((length(centers[1]), length(centers[2]), length(centers[3])))[cartesian_index]
                    push!(regular_cells, (cartesian_index, linear_index))
                end
            end
        end
    else
        error("Unsupported dimension: $dims")
    end

    return regular_cells
end

function identify!(MeshTag, mesh::CartesianMesh, body::Body)
    # Identify border cells
    border_cells = find_border(mesh)
    MeshTag.border_cells = border_cells
    
    # Identify cut cells
    cut_cells = find_cut(mesh, body)
    MeshTag.cut_cells = cut_cells

    # Identify regular cells
    regular_cells = find_regular(mesh)
    MeshTag.regular_cells = regular_cells

    return MeshTag
end
