function find_border(mesh)
    centers = mesh.centers
    dims = length(centers)
    border_cells = []

    if dims == 1
        for i in 1:length(centers[1])
            if i == 1 || i == length(centers[1])
                cartesian_index = CartesianIndex(i)
                linear_index = LinearIndices((length(centers[1])+1,))[cartesian_index]
                push!(border_cells, (cartesian_index, linear_index))
            end
        end
    elseif dims == 2
        for i in 1:length(centers[1])
            for j in 1:length(centers[2])
                if i == 1 || i == length(centers[1]) || j == 1 || j == length(centers[2])
                    cartesian_index = CartesianIndex(i, j)
                    linear_index = LinearIndices((length(centers[1])+1, length(centers[2])+1))[cartesian_index]
                    push!(border_cells, (cartesian_index, linear_index))
                end
            end
        end
    elseif dims == 3
        for i in 1:length(centers[1])
            for j in 1:length(centers[2])
                for k in 1:length(centers[3])
                    if i == 1 || i == length(centers[1]) || j == 1 || j == length(centers[2]) || k == 1 || k == length(centers[3])
                        cartesian_index = CartesianIndex(i, j, k)
                        linear_index = LinearIndices((length(centers[1])+1, length(centers[2])+1, length(centers[3])+1))[cartesian_index]
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

function identify!(mesh::CartesianMesh{N}, body::Body) where N
    # Identify border cells
    border_cells = find_border(mesh)
    mesh.tag.border_cells = border_cells
    
    # Identify cut cells
    cut_cells = find_cut(mesh, body)
    mesh.tag.cut_cells = cut_cells

    # Identify regular cells
    regular_cells = find_regular(mesh)
    mesh.tag.regular_cells = regular_cells

    return MeshTag
end

# Initialize temperature uniformly across the domain
function initialize_temperature_uniform!(T0ₒ::Vector{Float64}, T0ᵧ::Vector{Float64}, value::Float64)
    fill!(T0ₒ, value)
    fill!(T0ᵧ, value)
end

# Initialize temperature within a square region centered at 'center' with half-width 'half_width'
function initialize_temperature_square!(T0ₒ::Vector{Float64}, T0ᵧ::Vector{Float64}, x_coords::Vector{Float64}, y_coords::Vector{Float64}, center::Tuple{Float64, Float64}, half_width::Float64, value::Float64)
    nx = length(x_coords) - 1
    ny = length(y_coords) - 1
    for j in 1:(ny + 1)
        for i in 1:(nx + 1)
            idx = i + (j - 1) * (nx + 1)
            x = x_coords[i]
            y = y_coords[j]
            if abs(x - center[1]) <= half_width && abs(y - center[2]) <= half_width
                T0ₒ[idx] = value
                T0ᵧ[idx] = value
            end
        end
    end
end

# Initialize temperature within a circular region centered at 'center' with radius 'radius'
function initialize_temperature_circle!(T0ₒ::Vector{Float64}, T0ᵧ::Vector{Float64}, x_coords::Vector{Float64}, y_coords::Vector{Float64}, center::Tuple{Float64, Float64}, radius::Float64, value::Float64)
    nx = length(x_coords) - 1
    ny = length(y_coords) - 1
    for j in 1:(ny + 1)
        for i in 1:(nx + 1)
            idx = i + (j - 1) * (nx + 1)
            x = x_coords[i]
            y = y_coords[j]
            distance = sqrt((x - center[1])^2 + (y - center[2])^2)
            if distance <= radius
                T0ₒ[idx] = value
                T0ᵧ[idx] = value
            end
        end
    end
end

# Initialize temperature using a custom function 'func(x, y)'
function initialize_temperature_function!(T0ₒ::Vector{Float64}, T0ᵧ::Vector{Float64}, x_coords::Vector{Float64}, y_coords::Vector{Float64}, func::Function)
    nx = length(x_coords) - 1
    ny = length(y_coords) - 1
    for j in 1:(ny + 1)
        for i in 1:(nx + 1)
            idx = i + (j - 1) * (nx + 1)
            x = x_coords[i]
            y = y_coords[j]
            T_value = func(x, y)
            T0ₒ[idx] = T_value
            T0ᵧ[idx] = T_value
        end
    end
end