function plot_solution(solver, mesh::CartesianMesh, body::Body) # Déterminer la dimension du maillage 
    dims = length(mesh.centers)
    # Évaluer la fonction distance signée (SDF) sur les nœuds du maillage
    Z_sdf = if dims == 1
        [body.sdf(xi, 0.0, 0.0) for xi in mesh.nodes[1]]
    elseif dims == 2
        [body.sdf(xi, yi, 0.0) for yi in mesh.nodes[2], xi in mesh.nodes[1]]
    elseif dims == 3
        # Pour 3D, tracer une tranche centrale en z
        z_mid = div(length(mesh.nodes[3]), 2)
        [body.sdf(xi, yi, mesh.nodes[3][z_mid]) for yi in mesh.nodes[2], xi in mesh.nodes[1]]
    else
        error("Dimension non supportée: $dims")
    end

    # Déterminer le type de problème
    is_steady = isa(solver, DiffusionSteadyMono) || isa(solver, DiffusionSteadyDiph)
    is_monophasic = isa(solver, DiffusionSteadyMono) || isa(solver, DiffusionUnsteadyMono)

    # Tracer selon la dimension et le type de problème
    if dims == 1
        # Tracé en 1D
        if is_steady
            if is_monophasic
                u = solver.x
                x = mesh.centers[1]
                fig = Figure()
                ax = Axis(fig[1, 1], title="Monophasic Steady Solution", xlabel="x", ylabel="u")
                lines!(ax, x, u, color=:blue)
                # Tracer le zéro de SDF
                zero_indices = findall(u .== 0)
                scatter!(ax, x[zero_indices], zeros(length(zero_indices)), color=:red, label="SDF=0")
                legend!(ax)
                display(fig)
            elseif isa(solver, DiffusionSteadyDiph)
                u = solver.x
                u1 = u[1:length(u) ÷ 2]
                u2 = u[length(u) ÷ 2 + 1:end]
                x = mesh.centers[1]
                fig = Figure()
                ax = Axis(fig[1, 1], title="Diphasic Steady Solutions", xlabel="x", ylabel="u")
                lines!(ax, x, u1, color=:blue, label="Phase 1")
                lines!(ax, x, u2, color=:green, label="Phase 2")
                # Tracer le zéro de SDF
                zero_indices = findall(u1 .== 0)
                scatter!(ax, x[zero_indices], zeros(length(zero_indices)), color=:red, label="SDF=0")
                legend!(ax)
                display(fig)
            end
        else
            # Tracé unsteady en 1D
            if is_monophasic
                states = solver.states
                fig = Figure()
                ax = Axis(fig[1, 1], title="Monophasic Unsteady Solutions", xlabel="x", ylabel="u")
                for state in states
                    lines!(ax, mesh.centers[1], state, color=:blue, alpha=0.3)
                end
                # Tracer le zéro de SDF pour la première solution
                zero_indices = findall(states[1] .== 0)
                scatter!(ax, mesh.centers[1][zero_indices], zeros(length(zero_indices)), color=:red, label="SDF=0")
                legend!(ax)
                display(fig)
            elseif isa(solver, DiffusionUnsteadyDiph)
                states1 = [state[1] for state in solver.states]  # Phase 1
                states2 = [state[2] for state in solver.states]  # Phase 2
                fig = Figure()
                ax1 = Axis(fig[1, 1], title="Diphasic Unsteady - Phase 1", xlabel="x", ylabel="u1")
                ax2 = Axis(fig[1, 2], title="Diphasic Unsteady - Phase 2", xlabel="x", ylabel="u2")
                for (u1, u2) in zip(states1, states2)
                    lines!(ax1, mesh.centers[1], u1, color=:blue, alpha=0.3)
                    lines!(ax2, mesh.centers[1], u2, color=:green, alpha=0.3)
                end
                # Tracer le zéro de SDF pour la première solution
                zero_indices = findall(states1[1] .== 0)
                scatter!(ax1, mesh.centers[1][zero_indices], zeros(length(zero_indices)), color=:red, label="SDF=0")
                scatter!(ax2, mesh.centers[1][zero_indices], zeros(length(zero_indices)), color=:red, label="SDF=0")
                legend!(ax1)
                legend!(ax2)
                display(fig)
            end
        end

    elseif dims == 2
        # Tracé en 2D
        if is_steady
            fig = Figure(resolution=(800, 600))
            if is_monophasic
                u = solver.x
                reshaped_u = reshape(u, (length(mesh.centers[1]), length(mesh.centers[2])) )'
                ax = Axis(fig[1, 1], title="Monophasic Steady Diffusion", xlabel="x", ylabel="y")
                heatmap!(ax, mesh.centers[1], mesh.centers[2], reshaped_u, colormap=:viridis)
                contour!(ax, mesh.nodes[1], mesh.nodes[2], Z_sdf, levels=[0.0], color=:red, linewidth=2, label="SDF=0")
                Colorbar(fig[1, 2], ax.plots[1], label="Temperature")
            elseif isa(solver, DiffusionSteadyDiph)
                u1 = solver.x1
                u2 = solver.x2
                reshaped_u1 = reshape(u1, (length(mesh.centers[1]), length(mesh.centers[2])) )'
                reshaped_u2 = reshape(u2, (length(mesh.centers[1]), length(mesh.centers[2])) )'
                
                ax1 = Axis(fig[1, 1], title="Diphasic Steady - Phase 1", xlabel="x", ylabel="y")
                heatmap!(ax1, mesh.centers[1], mesh.centers[2], reshaped_u1, colormap=:viridis)
                contour!(ax1, mesh.nodes[1], mesh.nodes[2], Z_sdf, levels=[0.0], color=:black, linewidth=2, label="SDF=0")
                Colorbar(fig[1, 2], ax1.plots[1], label="Phase 1 Temperature")
                
                ax2 = Axis(fig[1, 3], title="Diphasic Steady - Phase 2", xlabel="x", ylabel="y")
                heatmap!(ax2, mesh.centers[1], mesh.centers[2], reshaped_u2, colormap=:viridis)
                contour!(ax2, mesh.nodes[1], mesh.nodes[2], Z_sdf, levels=[0.0], color=:black, linewidth=2, label="SDF=0")
                Colorbar(fig[1, 4], ax2.plots[1], label="Phase 2 Temperature")
            end
            display(fig)
        else
            # Tracé unsteady en 2D
            if is_monophasic
                states = solver.states
                reshaped_states = [reshape(state, (length(mesh.centers[1]), length(mesh.centers[2])) )' for state in states]
                min_val = minimum([minimum(u) for u in reshaped_states])
                max_val = maximum([maximum(u) for u in reshaped_states])
                
                fig = Figure(resolution=(800, 600))
                ax = Axis(fig[1, 1], title="Monophasic Unsteady Diffusion", xlabel="x", ylabel="y")
                hm = heatmap!(ax, mesh.centers[1], mesh.centers[2], reshaped_states[1], colormap=:viridis, colorrange=(min_val, max_val))
                contour!(ax, mesh.nodes[1], mesh.nodes[2], Z_sdf, levels=[0.0], color=:red, linewidth=2)
                Colorbar(fig[1, 2], hm, label="Temperature")
                
                # Animation
                record(fig, "heat_MonoUnsteady.mp4", 1:length(states); framerate=10) do frame
                    hm[1] = reshaped_states[frame]
                end
                display(fig)
                
            elseif isa(solver, DiffusionUnsteadyDiph)
                states = solver.states
                reshaped_u1 = [reshape(state[1], (length(mesh.centers[1]), length(mesh.centers[2])) )' for state in states]
                reshaped_u2 = [reshape(state[2], (length(mesh.centers[1]), length(mesh.centers[2])) )' for state in states]
                min_val = minimum([minimum(u) for u in reshaped_u1])
                max_val = maximum([maximum(u) for u in reshaped_u1])
                
                fig = Figure(resolution=(1600, 600))
                ax1 = Axis(fig[1, 1], title="Diphasic Unsteady - Phase 1", xlabel="x", ylabel="y")
                hm1 = heatmap!(ax1, mesh.centers[1], mesh.centers[2], reshaped_u1[1], colormap=:viridis, colorrange=(min_val, max_val))
                contour!(ax1, mesh.nodes[1], mesh.nodes[2], Z_sdf, levels=[0.0], color=:black, linewidth=2, label="SDF=0")
                Colorbar(fig[1, 2], hm1, label="Phase 1 Temperature")
                
                ax2 = Axis(fig[1, 3], title="Diphasic Unsteady - Phase 2", xlabel="x", ylabel="y")
                hm2 = heatmap!(ax2, mesh.centers[1], mesh.centers[2], reshaped_u2[1], colormap=:viridis, colorrange=(min_val, max_val))
                contour!(ax2, mesh.nodes[1], mesh.nodes[2], Z_sdf, levels=[0.0], color=:black, linewidth=2, label="SDF=0")
                Colorbar(fig[1, 4], hm2, label="Phase 2 Temperature")
                
                # Animation
                record(fig, "heat_DiffUnsteadyDiph.mp4", 1:length(states); framerate=10) do frame
                    hm1[1] = reshaped_u1[frame]
                    hm2[1] = reshaped_u2[frame]
                end
                display(fig)
            end
        end
        
    elseif dims == 3
        # Tracé en 3D (exemple avec une tranche centrale en z)
        if is_steady
            fig = Figure(resolution=(800, 600))
            if is_monophasic
                u = solver.x
                nx, ny, nz = length(mesh.centers[1]), length(mesh.centers[2]), length(mesh.centers[3])
                reshaped_u = reshape(u, (nx, ny, nz))
                z_mid = div(nz, 2)
                ax = Axis3(fig[1, 1], title="Monophasic Steady Diffusion (Slice)", xlabel="x", ylabel="y", zlabel="z")
                heatmap!(ax, mesh.centers[1], mesh.centers[2], reshaped_u[:, :, z_mid], colormap=:viridis)
                contour!(ax, mesh.nodes[1], mesh.nodes[2], mesh.nodes[3][z_mid], reshaped_u[:, :, z_mid], levels=[0.0], color=:red, linewidth=2)
                display(fig)
            elseif isa(solver, DiffusionSteadyDiph)
                u1 = solver.x1
                u2 = solver.x2
                nx, ny, nz = length(mesh.centers[1]), length(mesh.centers[2]), length(mesh.centers[3])
                reshaped_u1 = reshape(u1, (nx, ny, nz))
                reshaped_u2 = reshape(u2, (nx, ny, nz))
                z_mid = div(nz, 2)
                
                ax1 = Axis3(fig[1, 1], title="Diphasic Steady - Phase 1 (Slice)", xlabel="x", ylabel="y", zlabel="z")
                heatmap!(ax1, mesh.centers[1], mesh.centers[2], reshaped_u1[:, :, z_mid], colormap=:viridis)
                contour!(ax1, mesh.nodes[1], mesh.nodes[2], mesh.nodes[3][z_mid], reshaped_u1[:, :, z_mid], levels=[0.0], color=:black, linewidth=2)
                Colorbar(fig[1, 2], ax1.plots[1], label="Phase 1 Temperature")
                
                ax2 = Axis3(fig[1, 3], title="Diphasic Steady - Phase 2 (Slice)", xlabel="x", ylabel="y", zlabel="z")
                heatmap!(ax2, mesh.centers[1], mesh.centers[2], reshaped_u2[:, :, z_mid], colormap=:viridis)
                contour!(ax2, mesh.nodes[1], mesh.nodes[2], mesh.nodes[3][z_mid], reshaped_u2[:, :, z_mid], levels=[0.0], color=:black, linewidth=2)
                Colorbar(fig[1, 4], ax2.plots[1], label="Phase 2 Temperature")
                
                display(fig)
            end
        else
            # Tracé unsteady en 3D
            if is_monophasic
                states = solver.states
                nx, ny, nz = length(mesh.centers[1]), length(mesh.centers[2]), length(mesh.centers[3])
                reshaped_states = [reshape(state, (nx, ny, nz)) for state in states]
                z_mid = div(nz, 2)
                min_val = minimum([minimum(reshape(state[:, :, z_mid], :) ) for state in reshaped_states])
                max_val = maximum([maximum(reshape(state[:, :, z_mid], :) ) for state in reshaped_states])
                
                fig = Figure(resolution=(800, 600))
                ax = Axis3(fig[1, 1], title="Monophasic Unsteady Diffusion (Slice)", xlabel="x", ylabel="y", zlabel="z")
                hm = heatmap!(ax, mesh.centers[1], mesh.centers[2], reshaped_states[1][:, :, z_mid], colormap=:viridis, colorrange=(min_val, max_val))
                contour!(ax, mesh.nodes[1], mesh.nodes[2], mesh.nodes[3][z_mid], reshaped_states[1][:, :, z_mid], levels=[0.0], color=:red, linewidth=2)
                Colorbar(fig[1, 2], hm, label="Temperature")
                
                # Animation
                record(fig, "heat_MonoUnsteady3D.mp4", 1:length(states); framerate=10) do frame
                    hm[1] = reshaped_states[frame][:, :, z_mid]
                end
                display(fig)
                
            elseif isa(solver, DiffusionUnsteadyDiph)
                states = solver.states
                nx, ny, nz = length(mesh.centers[1]), length(mesh.centers[2]), length(mesh.centers[3])
                reshaped_u1 = [reshape(state[1], (nx, ny, nz)) for state in states]
                reshaped_u2 = [reshape(state[2], (nx, ny, nz)) for state in states]
                z_mid = div(nz, 2)
                min_val = minimum([minimum(reshape(u1[:, :, z_mid], :)) for u1 in reshaped_u1])
                max_val = maximum([maximum(reshape(u1[:, :, z_mid], :)) for u1 in reshaped_u1])
                
                fig = Figure(resolution=(1600, 600))
                ax1 = Axis3(fig[1, 1], title="Diphasic Unsteady - Phase 1 (Slice)", xlabel="x", ylabel="y", zlabel="z")
                hm1 = heatmap!(ax1, mesh.centers[1], mesh.centers[2], reshaped_u1[1][:, :, z_mid], colormap=:viridis, colorrange=(min_val, max_val))
                contour!(ax1, mesh.nodes[1], mesh.nodes[2], mesh.nodes[3][z_mid], reshaped_u1[1][:, :, z_mid], levels=[0.0], color=:black, linewidth=2)
                Colorbar(fig[1, 2], hm1, label="Phase 1 Temperature")
                
                ax2 = Axis3(fig[1, 3], title="Diphasic Unsteady - Phase 2 (Slice)", xlabel="x", ylabel="y", zlabel="z")
                hm2 = heatmap!(ax2, mesh.centers[1], mesh.centers[2], reshaped_u2[1][:, :, z_mid], colormap=:viridis, colorrange=(min_val, max_val))
                contour!(ax2, mesh.nodes[1], mesh.nodes[2], mesh.nodes[3][z_mid], reshaped_u2[1][:, :, z_mid], levels=[0.0], color=:black, linewidth=2)
                Colorbar(fig[1, 4], hm2, label="Phase 2 Temperature")
                
                # Animation
                record(fig, "heat_DiffUnsteadyDiph3D.mp4", 1:length(states); framerate=10) do frame
                    hm1[1] = reshaped_u1[frame][:, :, z_mid]
                    hm2[1] = reshaped_u2[frame][:, :, z_mid]
                end
                display(fig)
            end
        end
    else
        error("Dimension non supportée: $dims")
    end
