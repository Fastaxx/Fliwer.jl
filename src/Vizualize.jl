function plot_solution(solver, mesh::CartesianMesh, body::Body) # Déterminer la dimension du maillage 
    dims = length(mesh.centers)
    
    # Évaluer la fonction distance signée (SDF) sur les nœuds du maillage
    Z_sdf = if dims == 1
        [body.sdf(xi) for xi in mesh.nodes[1]]
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
    is_steady = solver.time_type == Steady # Problème stationnaire
    is_monophasic = solver.phase_type == Monophasic # Problème monophasique

    # Tracer selon la dimension et le type de problème
    if dims == 1
        # Tracé en 1D
        if is_steady
            if is_monophasic # Monophasic
                uₒ = solver.x[1:length(solver.x) ÷ 2]
                uᵧ = solver.x[length(solver.x) ÷ 2 + 1:end]
                x = mesh.centers[1]
                fig = Figure()
                ax = Axis(fig[1, 1], title="Monophasic Steady Solution", xlabel="x", ylabel="u")
                lines!(ax, uₒ, color=:blue, label="Bulk")
                lines!(ax, uᵧ, color=:green, label="Interface")
                axislegend(ax)
                display(fig)
            else # Diphasic
                u1ₒ = solver.x[1:length(solver.x) ÷ 4]
                u1ᵧ = solver.x[length(solver.x) ÷ 4 + 1:2*length(solver.x) ÷ 4]
                u2ₒ = solver.x[2*length(solver.x) ÷ 4 + 1:3*length(solver.x) ÷ 4]
                u2ᵧ = solver.x[3*length(solver.x) ÷ 4:end]
                x = mesh.centers[1]
                fig = Figure()
                ax = Axis(fig[1, 1], title="Diphasic Steady Solutions", xlabel="x", ylabel="u")
                lines!(ax, u1ₒ, color=:blue, label="Phase 1 - Bulk")
                lines!(ax, u1ᵧ, color=:green, label="Phase 1 - Interface")
                lines!(ax, u2ₒ, color=:green, label="Phase 2 - Bulk")
                lines!(ax, u2ᵧ, color=:blue, label="Phase 2 - Interface")
                axislegend(ax)
                display(fig)
            end
        else
            # Tracé unsteady en 1D
            if is_monophasic # Monophasic
                states = solver.states
                fig = Figure()
                ax = Axis(fig[1, 1], title="Monophasic Unsteady Solutions", xlabel="x", ylabel="u")
                for state in states
                    lines!(ax, mesh.centers[1], state[1:length(state) ÷ 2], color=:blue, alpha=0.3, label="Bulk")
                    lines!(ax, mesh.centers[1], state[length(state) ÷ 2 + 1:end], color=:green, alpha=0.3, label="Interface")
                end
                axislegend(ax)
                display(fig)
            else # Diphasic
                states1ₒ = [state[1:length(state) ÷ 4] for state in solver.states]  # Phase 1 - Bulk
                states1ᵧ = [state[length(state) ÷ 4 + 1:2*length(state) ÷ 4] for state in solver.states]  # Phase 1 - Interface
                states2ₒ = [state[2*length(state) ÷ 4 + 1:3*length(state) ÷ 4] for state in solver.states]  # Phase 2 - Bulk
                states2ᵧ = [state[3*length(state) ÷ 4:end] for state in solver.states]  # Phase 2 - Interface
                fig = Figure()
                ax1 = Axis(fig[1, 1], title="Diphasic Unsteady - Phase 1", xlabel="x", ylabel="u1")
                ax2 = Axis(fig[1, 2], title="Diphasic Unsteady - Phase 2", xlabel="x", ylabel="u2")
                for (state1ₒ, state1ᵧ, state2ₒ, state2ᵧ) in zip(states1ₒ, states1ᵧ, states2ₒ, states2ᵧ)
                    lines!(ax1, mesh.centers[1], state1ₒ, color=:blue, alpha=0.3, label="Bulk")
                    lines!(ax1, mesh.centers[1], state1ᵧ, color=:green, alpha=0.3, label="Interface")
                    lines!(ax2, mesh.centers[1], state2ₒ, color=:blue, alpha=0.3, label="Bulk")
                    lines!(ax2, mesh.centers[1], state2ᵧ, color=:green, alpha=0.3, label="Interface")
                end
                axislegend(ax1)
                axislegend(ax2)
                display(fig)
            end
        end
    elseif dims == 2
        # Tracé en 2D
        if is_steady
            fig = Figure(size=(800, 600))
            if is_monophasic # Monophasic
                uₒ = solver.x[1:length(solver.x) ÷ 2]
                uᵧ = solver.x[length(solver.x) ÷ 2 + 1:end]
                reshaped_uₒ = reshape(uₒ, (length(mesh.centers[1])+1, length(mesh.centers[2])+1) )'
                reshaped_uᵧ = reshape(uᵧ, (length(mesh.centers[1])+1, length(mesh.centers[2])+1) )'
                ax1 = Axis(fig[1, 1], title="Monophasic Steady Solution - Bulk", xlabel="x", ylabel="y", aspect = DataAspect())
                ax2 = Axis(fig[1, 3], title="Monophasic Steady Solution - Interface", xlabel="x", ylabel="y", aspect = DataAspect())
                hm1 = heatmap!(ax1, mesh.centers[1], mesh.centers[2], reshaped_uₒ, colormap=:viridis)
                hm2 = heatmap!(ax2, mesh.centers[1], mesh.centers[2], reshaped_uᵧ, colormap=:viridis)
                contour!(ax1, mesh.nodes[1], mesh.nodes[2], Z_sdf, levels=[0.0], color=:red, linewidth=2, label="SDF=0")
                contour!(ax2, mesh.nodes[1], mesh.nodes[2], Z_sdf, levels=[0.0], color=:red, linewidth=2, label="SDF=0")
                Colorbar(fig[1, 2], hm1, label="Bulk Temperature")
                Colorbar(fig[1, 4], hm2, label="Interface Temperature")
            else # Diphasic
                u1ₒ = solver.x[1:length(solver.x) ÷ 4]
                u1ᵧ = solver.x[length(solver.x) ÷ 4 + 1:2*length(solver.x) ÷ 4]
                u2ₒ = solver.x[2*length(solver.x) ÷ 4 + 1:3*length(solver.x) ÷ 4]
                u2ᵧ = solver.x[3*length(solver.x) ÷ 4+1:end]
                reshaped_u1ₒ = reshape(u1ₒ, (length(mesh.centers[1])+1, length(mesh.centers[2])+1) )'
                reshaped_u1ᵧ = reshape(u1ᵧ, (length(mesh.centers[1])+1, length(mesh.centers[2])+1) )'
                reshaped_u2ₒ = reshape(u2ₒ, (length(mesh.centers[1])+1, length(mesh.centers[2])+1) )'
                reshaped_u2ᵧ = reshape(u2ᵧ, (length(mesh.centers[1])+1, length(mesh.centers[2])+1) )'
                
                ax1 = Axis(fig[1, 1], title="Diphasic Steady - Phase 1 - Bulk", xlabel="x", ylabel="y", aspect = DataAspect())
                hm1 = heatmap!(ax1, mesh.centers[1], mesh.centers[2], reshaped_u1ₒ, colormap=:viridis)
                cb1 = Colorbar(fig[1, 2], hm1, label="Phase 1 Bulk Temperature")
                contour!(ax1, mesh.nodes[1], mesh.nodes[2], Z_sdf, levels=[0.0], color=:red, linewidth=2, label="SDF=0")

                ax2 = Axis(fig[1, 3], title="Diphasic Steady - Phase 1 - Interface", xlabel="x", ylabel="y", aspect = DataAspect())
                hm2 = heatmap!(ax2, mesh.centers[1], mesh.centers[2], reshaped_u1ᵧ, colormap=:viridis)
                cb2 = Colorbar(fig[1, 4], hm2, label="Phase 1 Interface Temperature")
                contour!(ax2, mesh.nodes[1], mesh.nodes[2], Z_sdf, levels=[0.0], color=:red, linewidth=2, label="SDF=0")

                ax3 = Axis(fig[2, 1], title="Diphasic Steady - Phase 2 - Bulk", xlabel="x", ylabel="y", aspect = DataAspect())
                hm3 = heatmap!(ax3, mesh.centers[1], mesh.centers[2], reshaped_u2ₒ, colormap=:viridis)
                cb3 = Colorbar(fig[2, 2], hm3, label="Phase 2 Bulk Temperature")
                contour!(ax3, mesh.nodes[1], mesh.nodes[2], Z_sdf, levels=[0.0], color=:red, linewidth=2, label="SDF=0")

                ax4 = Axis(fig[2, 3], title="Diphasic Steady - Phase 2 - Interface", xlabel="x", ylabel="y", aspect = DataAspect())
                hm4 = heatmap!(ax4, mesh.centers[1], mesh.centers[2], reshaped_u2ᵧ, colormap=:viridis)
                cb4 = Colorbar(fig[2, 4], hm4, label="Phase 2 Interface Temperature")
                contour!(ax4, mesh.nodes[1], mesh.nodes[2], Z_sdf, levels=[0.0], color=:red, linewidth=2, label="SDF=0")

            end
            display(fig)
        else
            # Tracé unsteady en 2D
            if is_monophasic
                uₒ = solver.states[1][1:length(solver.states[1]) ÷ 2]
                uᵧ = solver.states[1][length(solver.states[1]) ÷ 2 + 1:end]

                reshaped_uₒ = reshape(uₒ, (length(mesh.centers[1])+1, length(mesh.centers[2])+1) )'
                reshaped_uᵧ = reshape(uᵧ, (length(mesh.centers[1])+1, length(mesh.centers[2])+1) )'

                fig = Figure(size=(800, 600))

                ax1 = Axis(fig[1, 1], title="Monophasic Unsteady Diffusion - Bulk", xlabel="x", ylabel="y", aspect = DataAspect())
                hm1 = heatmap!(ax1, mesh.centers[1], mesh.centers[2], reshaped_uₒ, colormap=:viridis)
                contour!(ax1, mesh.nodes[1], mesh.nodes[2], Z_sdf, levels=[0.0], color=:red, linewidth=2, label="SDF=0")
                Colorbar(fig[1, 2], hm1, label="Bulk Temperature")

                ax2 = Axis(fig[1, 3], title="Monophasic Unsteady Diffusion - Interface", xlabel="x", ylabel="y", aspect = DataAspect())
                hm2 = heatmap!(ax2, mesh.centers[1], mesh.centers[2], reshaped_uᵧ, colormap=:viridis)
                contour!(ax2, mesh.nodes[1], mesh.nodes[2], Z_sdf, levels=[0.0], color=:red, linewidth=2, label="SDF=0")
                Colorbar(fig[1, 4], hm2, label="Interface Temperature")

                display(fig)

                # Animation
                fig = Figure(size=(800, 400))
                ax = Axis(fig[1, 1], title="Diffusion Unsteady Monophasic", xlabel="x", ylabel="y", aspect = DataAspect())
                xlims!(ax, mesh.x0[1], length(mesh.centers[1]))
                ylims!(ax, mesh.x0[2], length(mesh.centers[2]))

                min_val = minimum([minimum(reshape(state[1:length(state) ÷ 2], (length(mesh.centers[1])+1, length(mesh.centers[2])+1))') for state in solver.states])
                max_val = maximum([maximum(reshape(state[1:length(state) ÷ 2], (length(mesh.centers[1])+1, length(mesh.centers[2])+1))') for state in solver.states])

                hm = heatmap!(ax, reshape(solver.states[1][1:length(solver.states[1]) ÷ 2], (length(mesh.centers[1])+1, length(mesh.centers[2])+1))', colormap=:viridis, colorrange=(min_val, max_val))
                Colorbar(fig[1, 2], hm, label="Temperature")

                update_hm(frame) = reshape(solver.states[frame][1:length(solver.states[frame]) ÷ 2], (length(mesh.centers[1])+1, length(mesh.centers[2])+1))'

                record(fig, "heat_MonoUnsteady.mp4", 1:length(solver.states); framerate=10) do frame
                    hm[1] = update_hm(frame)
                end

                display(fig)
                
            else
                # Plot Last State
                u1ₒ = solver.states[1][1:length(solver.states[1]) ÷ 4]
                u1ᵧ = solver.states[1][length(solver.states[1]) ÷ 4 + 1:2*length(solver.states[1]) ÷ 4]
                u2ₒ = solver.states[1][2*length(solver.states[1]) ÷ 4 + 1:3*length(solver.states[1]) ÷ 4]
                u2ᵧ = solver.states[1][3*length(solver.states[1]) ÷ 4 + 1:end]

                reshaped_u1ₒ = reshape(u1ₒ, (length(mesh.centers[1])+1, length(mesh.centers[2])+1) )'
                reshaped_u1ᵧ = reshape(u1ᵧ, (length(mesh.centers[1])+1, length(mesh.centers[2])+1) )'
                reshaped_u2ₒ = reshape(u2ₒ, (length(mesh.centers[1])+1, length(mesh.centers[2])+1) )'
                reshaped_u2ᵧ = reshape(u2ᵧ, (length(mesh.centers[1])+1, length(mesh.centers[2])+1) )'

                fig = Figure()

                x, y = range(mesh.x0[1], stop=mesh.x0[1]+mesh.h[1][1]*length(mesh.h[1]), length=length(mesh.h[1])+1), range(mesh.x0[2], stop=mesh.x0[2]+mesh.h[2][1]*length(mesh.h[2]), length=length(mesh.h[2])+1)

                ax1 = Axis(fig[1, 1], title="Diphasic Unsteady - Phase 1 - Bulk", xlabel="x", ylabel="y")
                surface!(ax1, x, y, reshaped_u1ₒ, colormap=:viridis)

                ax2 = Axis(fig[1, 2], title="Diphasic Unsteady - Phase 1 - Interface", xlabel="x", ylabel="y")
                surface!(ax2, x, y, reshaped_u1ᵧ, colormap=:viridis)

                ax3 = Axis(fig[2, 1], title="Diphasic Unsteady - Phase 2 - Bulk", xlabel="x", ylabel="y")
                surface!(ax3, x, y, reshaped_u2ₒ, colormap=:viridis)

                ax4 = Axis(fig[2, 2], title="Diphasic Unsteady - Phase 2 - Interface", xlabel="x", ylabel="y")
                surface!(ax4, x, y, reshaped_u2ᵧ, colormap=:viridis)

                display(fig)

                # Animation
                fig = Figure(size=(800, 400))
                
                ax1 = Axis3(fig[1, 1], title="Diphasic Unsteady - Phase 1 - Bulk", xlabel="x", ylabel="y", zlabel="Temperature")
                s1 = surface!(ax1, reshaped_u1ₒ, colormap=:viridis)

                ax2 = Axis3(fig[1, 2], title="Diphasic Unsteady - Phase 1 - Interface", xlabel="x", ylabel="y", zlabel="Temperature")
                s2 = surface!(ax2, reshaped_u1ᵧ, colormap=:viridis)

                ax3 = Axis3(fig[2, 1], title="Diphasic Unsteady - Phase 2 - Bulk", xlabel="x", ylabel="y", zlabel="Temperature")
                s3 = surface!(ax3, reshaped_u2ₒ, colormap=:viridis)

                ax4 = Axis3(fig[2, 2], title="Diphasic Unsteady - Phase 2 - Interface", xlabel="x", ylabel="y", zlabel="Temperature")
                s4 = surface!(ax4, reshaped_u2ᵧ, colormap=:viridis)

                zlims!(ax1, 0, 1)
                zlims!(ax2, 0, 1)
                zlims!(ax3, 0, 1)
                zlims!(ax4, 0, 1)

                function update_surfaces!(frame)
                    s1[:z] = reshape(solver.states[frame][1:length(solver.states[frame]) ÷ 4], (length(mesh.centers[1])+1, length(mesh.centers[2])+1))'
                    #s2[:z] = reshape(solver.states[frame][length(solver.states[frame]) ÷ 4 + 1:2*length(solver.states[frame]) ÷ 4], length(mesh.centers[1])+1, length(mesh.centers[2])+1))'
                    s3[:z] = reshape(solver.states[frame][2*length(solver.states[frame]) ÷ 4 + 1:3*length(solver.states[frame]) ÷ 4], (length(mesh.centers[1])+1, length(mesh.centers[2])+1))'
                    #s4[:z] = reshape(solver.states[frame][3*length(solver.states[frame]) ÷ 4 + 1:end], length(mesh.centers[1])+1, length(mesh.centers[2])+1))'
                    println("Frame $frame")
                end

                record(fig, "heat_DiphUnsteady.mp4", 1:length(solver.states); framerate=10) do frame
                    update_surfaces!(frame)
                end

                display(fig)
            end
        end
        
    elseif dims == 3
        # Tracé en 3D (exemple avec une tranche centrale en z)
        if is_steady
            fig = Figure(size=(800, 600))
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
                
                fig = Figure(size=(800, 600))
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
                
                fig = Figure(size=(1600, 600))
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

end

