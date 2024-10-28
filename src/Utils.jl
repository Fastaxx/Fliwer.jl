mutable struct MeshTag
    border_cells::Array{Int64, 1}
    cut_cells::Array{Int64, 1}
    regular_cells::Array{Int64, 1}
    intersection_points::Array{Int64, 1}
end

function find_border(mesh)
    centers = mesh.centers
    dims = length(centers)
    border_cells = []

    if dims == 1
        for i in 1:length(centers[1])
            if i == 1 || i == length(centers[1])
                cartesian_index = CartesianIndex(i)
                linear_index = i
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
        for i in 1:length(nodes[1])
            sdf[i] = body.sdf(nodes[1][i])
        end
    elseif dim == 2
        for j in 1:length(nodes[2])
            for i in 1:length(nodes[1])
                sdf[LinearIndices((length(nodes[1]), length(nodes[2])))[CartesianIndex(i, j)]] = body.sdf(nodes[1][i], nodes[2][j])
            end
        end
    elseif dim == 3
        for k in 1:length(nodes[3])
            for j in 1:length(nodes[2])
                for i in 1:length(nodes[1])
                    sdf[LinearIndices((length(nodes[1]), length(nodes[2]), length(nodes[3])))[CartesianIndex(i, j, k)]] = body.sdf(nodes[1][i], nodes[2][j], nodes[3][k])
                end
            end
        end
    else
        error("Unsupported dimension: $dim")
    end

    return sdf
end

function find_cut(mesh::CartesianMesh, body::Body)
    # Evaluate the signed distance function on the mesh
    sdf = eval_sdf(mesh, body)
    cut_cells = []
    dims = length(mesh.centers)
    # If the signed distance function changes sign in a cell, it is a cut cell
    if dims == 1
        for i in 1:length(sdf)-1
            if sign(sdf[i]) != sign(sdf[i+1])
                cartesian_index = CartesianIndex(i)
                linear_index = LinearIndices((length(mesh.centers[1]),))[cartesian_index]
                push!(cut_cells, (cartesian_index, linear_index))
            end
        end
    elseif dims == 2
        for j in 1:length(mesh.centers[2])
            for i in 1:length(mesh.centers[1])
                if sign(sdf[LinearIndices((length(mesh.centers[1]), length(mesh.centers[2])))[CartesianIndex(i, j)]] * sdf[LinearIndices((length(mesh.centers[1]), length(mesh.centers[2])))[CartesianIndex(i+1, j)]]) < 0 || 
                   sign(sdf[LinearIndices((length(mesh.centers[1]), length(mesh.centers[2])))[CartesianIndex(i, j)]] * sdf[LinearIndices((length(mesh.centers[1]), length(mesh.centers[2])))[CartesianIndex(i, j+1)]]) < 0
                    cartesian_index = CartesianIndex(i, j)
                    linear_index = LinearIndices((length(mesh.centers[1]), length(mesh.centers[2])))[cartesian_index]
                    push!(cut_cells, (cartesian_index, linear_index))
                end
            end
        end
    elseif dims == 3
        for k in 1:length(mesh.centers[3])-1
            for j in 1:length(mesh.centers[2])-1
                for i in 1:length(mesh.centers[1])-1
                    if sign(sdf[LinearIndices((length(mesh.centers[1]), length(mesh.centers[2]), length(mesh.centers[3])))[CartesianIndex(i, j, k)]] * sdf[LinearIndices((length(mesh.centers[1]), length(mesh.centers[2]), length(mesh.centers[3])))[CartesianIndex(i+1, j, k)]]) < 0 || 
                       sign(sdf[LinearIndices((length(mesh.centers[1]), length(mesh.centers[2]), length(mesh.centers[3])))[CartesianIndex(i, j, k)]] * sdf[LinearIndices((length(mesh.centers[1]), length(mesh.centers[2]), length(mesh.centers[3])))[CartesianIndex(i, j+1, k)]]) < 0 || 
                       sign(sdf[LinearIndices((length(mesh.centers[1]), length(mesh.centers[2]), length(mesh.centers[3])))[CartesianIndex(i, j, k)]] * sdf[LinearIndices((length(mesh.centers[1]), length(mesh.centers[2]), length(mesh.centers[3])))[CartesianIndex(i, j, k+1)]]) < 0
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


function identify!(MeshTag, mesh::CartesianMesh)
    # Identify border cells
    border_cells = find_border(mesh)
    MeshTag.border_cells = border_cells
    
    # Identify cut cells
    cut_cells = find_cut(mesh, body)
    MeshTag.cut_cells = cut_cells

    # Identify regular cells
    regular_cells = find_regular(mesh)
    MeshTag.regular_cells = regular_cells

    # Identify intersection points
    intersection_points = find_intersection(mesh)
    MeshTag.intersection_points = intersection_points

    return MeshTag
end
