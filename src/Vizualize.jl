function animate_solution(solver, mesh::CartesianMesh{1}, body::Body)
    # Déterminer le type de problème
    is_monophasic = solver.phase_type == Monophasic # Problème monophasique

    # Enregistrer l'animation selon le type de problème
    if is_monophasic
        # Récupérer les états
        states = solver.states

        # Créer une figure
        fig = Figure()

        # Créer un axe pour la figure
        ax = Axis(fig[1, 1], title="Monophasic Unsteady Diffusion", xlabel="x", ylabel="u")

        function update_ln(frame)
            # Récupérer l'état
            state = states[frame]

            # Tracer l'état
            lines!(ax, state[1:length(state) ÷ 2], color=:blue, alpha=0.3, label="Bulk")
            lines!(ax, state[length(state) ÷ 2 + 1:end], color=:green, alpha=0.3, label="Interface")
        end

        # Enregistrer l'animation
        record(fig, "heat_MonoUnsteady.mp4", 1:length(states); framerate=10) do frame
            update_ln(frame)
        end

        # Afficher la figure
        display(fig)
    else
        # Récupérer les états
        states = solver.states

        # Créer une figure
        fig = Figure()

        # Créer un axe pour la figure
        ax1 = Axis(fig[1, 1], title="Diphasic Unsteady - Phase 1", xlabel="x", ylabel="u1")
        ax2 = Axis(fig[2, 1], title="Diphasic Unsteady - Phase 2", xlabel="x", ylabel="u2")

        # Créer une fonction pour mettre à jour la figure
        function update_plot(frame)
            # Récupérer l'état
            state = states[frame]

            # Tracer l'état
            lines!(ax1, state[1:length(state) ÷ 4], color=:blue, alpha=0.3, label="Bulk")
            lines!(ax1, state[length(state) ÷ 4 + 1:2*length(state) ÷ 4], color=:green, alpha=0.3, label="Interface")
            lines!(ax2, state[2*length(state) ÷ 4 + 1:3*length(state) ÷ 4], color=:blue, alpha=0.3, label="Bulk")
            lines!(ax2, state[3*length(state) ÷ 4 + 1:end], color=:green, alpha=0.3, label="Interface")

        end

        # Enregistrer l'animation
        record(fig, "heat_DiphUnsteady.mp4", 1:length(states); framerate=10) do frame
            update_plot(frame)
        end

        # Afficher la figure
        display(fig)

    end
end

function animate_solution(solver, mesh::CartesianMesh{2}, body::Body)
    # Déterminer le type de problème
    is_monophasic = solver.phase_type == Monophasic # Problème monophasique

    # Enregistrer l'animation selon le type de problème
    if is_monophasic
        # Récupérer les états
        states = solver.states

        # Créer une figure
        fig = Figure()

        # Créer un axe pour la figure
        ax = Axis(fig[1, 1], title="Monophasic Unsteady", xlabel="x", ylabel="y", aspect=DataAspect())
        xlims!(ax, mesh.x0[1], length(mesh.centers[1]))
        ylims!(ax, mesh.x0[2], length(mesh.centers[2]))

        min_val = minimum([minimum(reshape(state[1:length(state) ÷ 2], (length(mesh.centers[1])+1, length(mesh.centers[2])+1))') for state in solver.states])
        max_val = maximum([maximum(reshape(state[1:length(state) ÷ 2], (length(mesh.centers[1])+1, length(mesh.centers[2])+1))') for state in solver.states])

        hm = heatmap!(ax, reshape(states[1][1:length(states[1]) ÷ 2], (length(mesh.centers[1])+1, length(mesh.centers[2])+1))', colormap=:viridis, colorrange=(min_val, max_val))
        Colorbar(fig[1, 2], hm, label="Temperature")

        update_hm(frame) = reshape(solver.states[frame][1:length(solver.states[frame]) ÷ 2], (length(mesh.centers[1])+1, length(mesh.centers[2])+1))'

        record(fig, "heat_MonoUnsteady.mp4", 1:length(solver.states); framerate=10) do frame
            hm[1] = update_hm(frame)
        end

        display(fig)
    else
        # Récupérer les états
        states = solver.states

        # Créer une figure
        fig = Figure(size=(800, 400))

        # Créer un axe pour la figure
        ax1 = Axis3(fig[1, 1], title="Diphasic Unsteady - Phase 1 - Bulk", xlabel="x", ylabel="y", zlabel="Temperature")
        s1 = surface!(ax1, reshape(states[1][1:length(states[1]) ÷ 4], (length(mesh.centers[1])+1, length(mesh.centers[2])+1))', colormap=:viridis)

        ax2 = Axis3(fig[1, 2], title="Diphasic Unsteady - Phase 1 - Interface", xlabel="x", ylabel="y", zlabel="Temperature")
        s2 = surface!(ax2, reshape(states[1][length(states[1]) ÷ 4 + 1:2*length(states[1]) ÷ 4], (length(mesh.centers[1])+1, length(mesh.centers[2])+1))', colormap=:viridis)

        ax3 = Axis3(fig[2, 1], title="Diphasic Unsteady - Phase 2 - Bulk", xlabel="x", ylabel="y", zlabel="Temperature")
        s3 = surface!(ax3, reshape(states[1][2*length(states[1]) ÷ 4 + 1:3*length(states[1]) ÷ 4], (length(mesh.centers[1])+1, length(mesh.centers[2])+1))', colormap=:viridis)

        ax4 = Axis3(fig[2, 2], title="Diphasic Unsteady - Phase 2 - Interface", xlabel="x", ylabel="y", zlabel="Temperature")
        s4 = surface!(ax4, reshape(states[1][3*length(states[1]) ÷ 4 + 1:end], (length(mesh.centers[1])+1, length(mesh.centers[2])+1))', colormap=:viridis)

        zlims!(ax1, 0, 1)
        zlims!(ax2, 0, 1)
        zlims!(ax3, 0, 1)
        zlims!(ax4, 0, 1)

        function update_surfaces!(frame)
            s1[:z] = reshape(solver.states[frame][1:length(solver.states[frame]) ÷ 4], (length(mesh.centers[1])+1, length(mesh.centers[2])+1))'
            s2[:z] = reshape(solver.states[frame][length(solver.states[frame]) ÷ 4 + 1:2*length(solver.states[frame]) ÷ 4], (length(mesh.centers[1])+1, length(mesh.centers[2])+1))'
            s3[:z] = reshape(solver.states[frame][2*length(solver.states[frame]) ÷ 4 + 1:3*length(solver.states[frame]) ÷ 4], (length(mesh.centers[1])+1, length(mesh.centers[2])+1))'
            s4[:z] = reshape(solver.states[frame][3*length(solver.states[frame]) ÷ 4 + 1:end], (length(mesh.centers[1])+1, length(mesh.centers[2])+1))'
            println("Frame $frame")
        end

        record(fig, "heat_DiphUnsteady.mp4", 1:length(solver.states); framerate=10) do frame
            update_surfaces!(frame)
        end

        display(fig)
    end
end

function animate_solution(solver, mesh::CartesianMesh{3}, body::Body)
    # Déterminer le type de problème
    is_monophasic = solver.phase_type == Monophasic # Problème monophasique

    if is_monophasic
        println("Animation Monophasic Unsteady Diffusion in 3D")
        # Initialisation de la figure et de la scène
        fig = Figure(size = (800, 600))
        ax = LScene(fig[1, 1], show_axis=false)

        # Définition des axes en fonction du maillage
        nx, ny, nz = length(mesh.centers[1]), length(mesh.centers[2]), length(mesh.centers[3])
        x = LinRange(mesh.x0[1], mesh.x0[1] + mesh.h[1][1]*nx, nx+1)
        y = LinRange(mesh.x0[2], mesh.x0[2] + mesh.h[2][1]*ny, ny+1)
        z = LinRange(mesh.x0[3], mesh.x0[3] + mesh.h[3][1]*nz, nz+1)
        z_index = Observable{Any}(1)
        # Création des sliders pour les plans de coupe
        sgrid = SliderGrid(
            fig[2, 1],
            (label = "Plan YZ - Axe X", range = 1:length(x)),
            (label = "Plan XZ - Axe Y", range = 1:length(y)),
            (label = "Plan XY - Axe Z", range = 1:length(z)),
            (label = "Frame", range = 1:length(solver.states)),
        )

        lo = sgrid.layout
        nc = ncols(lo)

        # Initialisation avec le premier état
        state_i = 1
        vol = solver.states[state_i][1:(end ÷ 2)]
        vol = reshape(vol, nx+1, ny+1, nz+1)
        plt = volumeslices!(ax, x, y, z, vol, colormap = :viridis)

        # Connexion des sliders aux méthodes de mise à jour des volumeslices
        sl_yz, sl_xz, sl_xy, sl_t = sgrid.sliders

        on(sl_yz.value) do v
            plt[:update_yz][](v)
        end

        on(sl_xz.value) do v
            plt[:update_xz][](v)
        end

        on(sl_xy.value) do v
            plt[:update_xy][](v)
        end

        on(sl_t.value) do v
            state_i = v
            vol = solver.states[state_i][1:(end ÷ 2)]
            vol = reshape(vol, nx+1, ny+1, nz+1)
            plt[:update_xy][](vol)
            plt[:update_xz][](vol)
            plt[:update_yz][](vol) 
        end

        connect!(z_index, sl_t.value)

        # Positionnement initial des sliders au centre
        set_close_to!(sl_yz, 0.5 * length(x))
        set_close_to!(sl_xz, 0.5 * length(y))
        set_close_to!(sl_xy, 0.5 * length(z))
        set_close_to!(sl_t, 1.0)

        # Ajout de toggles pour afficher/masquer les heatmaps
        hmaps = [plt[Symbol(:heatmap_, s)][] for s ∈ (:yz, :xz, :xy)]
        toggles = [Toggle(lo[i, nc + 1], active = true) for i ∈ 1:length(hmaps)]

        for (h, t) in zip(hmaps, toggles)
            connect!(h.visible, t.active)
        end
        
        display(fig)
    else 
        println("Diphasic Unsteady Diffusion in 3D")
    end

end



function plot_solution(solver, mesh::CartesianMesh{1}, body::Body, capacity::Capacity; state_i=1)
    # Déterminer le type de problème
    is_steady = solver.time_type == Steady # Problème stationnaire
    is_monophasic = solver.phase_type == Monophasic # Problème monophasique

    # Tracer selon la dimension et le type de problème
    if is_steady
        if is_monophasic # Monophasic
            uₒ = solver.x[1:length(solver.x) ÷ 2]
            uᵧ = solver.x[length(solver.x) ÷ 2 + 1:end]
            x = mesh.centers[1]

            # Désactiver les cellules désactivées
            cell_types = capacity.cell_types
            uₒ[cell_types .== 0] .= NaN
            #uᵧ[cell_types .== 0] .= NaN
            uᵧ[cell_types .== 1] .= NaN
            
            fig = Figure()
            ax = Axis(fig[1, 1], title="Monophasic Steady Solution", xlabel="x", ylabel="u")
            scatter!(ax, uₒ, color=:blue, label="Bulk")
            scatter!(ax, uᵧ, color=:green, label="Interface")
            axislegend(ax)
            display(fig)
        else # Diphasic
            u1ₒ = solver.x[1:length(solver.x) ÷ 4]
            u1ᵧ = solver.x[length(solver.x) ÷ 4 + 1:2*length(solver.x) ÷ 4]
            u2ₒ = solver.x[2*length(solver.x) ÷ 4 + 1:3*length(solver.x) ÷ 4]
            u2ᵧ = solver.x[3*length(solver.x) ÷ 4 + 1:end]

            # Désactiver les cellules désactivées
            cell_types = capacity.cell_types
            u1ₒ[cell_types .== 0] .= NaN
            u1ᵧ[cell_types .== 0] .= NaN
            u1ᵧ[cell_types .== 1] .= NaN
            u2ₒ[cell_types .== 1] .= NaN
            u2ᵧ[cell_types .== 1] .= NaN
            u2ᵧ[cell_types .== 0] .= NaN

            x = mesh.centers[1]
            fig = Figure()
            ax = Axis(fig[1, 1], title="Diphasic Steady Solutions", xlabel="x", ylabel="u")
            scatter!(ax, u1ₒ, color=:blue, label="Phase 1 - Bulk")
            scatter!(ax, u1ᵧ, color=:green, label="Phase 1 - Interface")
            scatter!(ax, u2ₒ, color=:red, label="Phase 2 - Bulk")
            scatter!(ax, u2ᵧ, color=:purple, label="Phase 2 - Interface")
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
                lines!(ax, state[1:length(state) ÷ 2], color=:blue, alpha=0.3, label="Bulk")
                lines!(ax, state[length(state) ÷ 2 + 1:end], color=:green, alpha=0.3, label="Interface")
            end
            #axislegend(ax)
            display(fig)
        else # Diphasic

            states1ₒ = solver.states[state_i][1:length(solver.states[state_i]) ÷ 4]  # Phase 1 - Bulk
            states1ᵧ = solver.states[state_i][length(solver.states[state_i]) ÷ 4 + 1:2*length(solver.states[state_i]) ÷ 4]  # Phase 1 - Interface
            states2ₒ = solver.states[state_i][2*length(solver.states[state_i]) ÷ 4 + 1:3*length(solver.states[state_i]) ÷ 4]  # Phase 2 - Bulk
            states2ᵧ = solver.states[state_i][3*length(solver.states[state_i]) ÷ 4 + 1:end]  # Phase 2 - Interface

            cell_types = capacity.cell_types
            states1ₒ[cell_types .== 0] .= NaN
            states1ᵧ[cell_types .== 0] .= NaN
            states1ᵧ[cell_types .== 1] .= NaN
            states2ₒ[cell_types .== 1] .= NaN
            states2ᵧ[cell_types .== 1] .= NaN
            states2ᵧ[cell_types .== 0] .= NaN

            fig = Figure(size=(800, 600))
            ax1 = Axis(fig[1, 1], title="Diphasic Unsteady - Phase 1", xlabel="x", ylabel="u1")
            ax2 = Axis(fig[2, 1], title="Diphasic Unsteady - Phase 2", xlabel="x", ylabel="u2")
            lines!(ax1, states1ₒ, color=:blue, label="Bulk")
            scatter!(ax1, states1ᵧ, color=:green, label="Interface")
            lines!(ax2, states2ₒ, color=:blue, label="Bulk")
            scatter!(ax2, states2ᵧ, color=:green, label="Interface")
            #axislegend(ax1)
            #axislegend(ax2)
            display(fig)
        end
    end
end

function plot_solution(solver, mesh::CartesianMesh{2}, body::Body, capacity::Capacity; state_i=1)
    Z_sdf = [body.sdf(xi, yi, 0.0) for yi in mesh.nodes[2], xi in mesh.nodes[1]]

    is_steady = solver.time_type == Steady
    is_monophasic = solver.phase_type == Monophasic

    # Tracé en 2D
    if is_steady
        fig = Figure(size=(800, 600))
        if is_monophasic # Monophasic
            cell_types = capacity.cell_types
            uₒ = solver.x[1:length(solver.x) ÷ 2]
            uᵧ = solver.x[length(solver.x) ÷ 2 + 1:end]

            # Désactiver les cellules désactivées
            uₒ[cell_types .== 0] .= NaN
            uᵧ[cell_types .== 0] .= NaN
            uᵧ[cell_types .== 1] .= NaN

            reshaped_uₒ = reshape(uₒ, (length(mesh.centers[1])+1, length(mesh.centers[2])+1) )'
            reshaped_uᵧ = reshape(uᵧ, (length(mesh.centers[1])+1, length(mesh.centers[2])+1) )'

            ax1 = Axis(fig[1, 1], title="Monophasic Steady Solution - Bulk", xlabel="x", ylabel="y", aspect = DataAspect())
            ax2 = Axis(fig[1, 3], title="Monophasic Steady Solution - Interface", xlabel="x", ylabel="y", aspect = DataAspect())

            hm1 = heatmap!(ax1, mesh.centers[1], mesh.centers[2], reshaped_uₒ, colormap=:viridis)
            hm2 = heatmap!(ax2, mesh.centers[1], mesh.centers[2], reshaped_uᵧ, colormap=:viridis)

            #contour!(ax1, mesh.nodes[1], mesh.nodes[2], Z_sdf, levels=[0.0], color=:red, linewidth=2, label="SDF=0")
            #contour!(ax2, mesh.nodes[1], mesh.nodes[2], Z_sdf, levels=[0.0], color=:red, linewidth=2, label="SDF=0")

            Colorbar(fig[1, 2], hm1, label="Bulk Temperature")
            Colorbar(fig[1, 4], hm2, label="Interface Temperature")
        else # Diphasic
            cell_types = capacity.cell_types
            u1ₒ = solver.x[1:length(solver.x) ÷ 4]
            u1ᵧ = solver.x[length(solver.x) ÷ 4 + 1:2*length(solver.x) ÷ 4]
            u2ₒ = solver.x[2*length(solver.x) ÷ 4 + 1:3*length(solver.x) ÷ 4]
            u2ᵧ = solver.x[3*length(solver.x) ÷ 4+1:end]

            # Désactiver les cellules désactivées
            u1ₒ[cell_types .== 0] .= NaN
            u1ᵧ[cell_types .== 0] .= NaN
            u1ᵧ[cell_types .== 1] .= NaN
            u2ₒ[cell_types .== 1] .= NaN
            u2ᵧ[cell_types .== 1] .= NaN
            u2ᵧ[cell_types .== 0] .= NaN

            reshaped_u1ₒ = reshape(u1ₒ, (length(mesh.centers[1])+1, length(mesh.centers[2])+1) )'
            reshaped_u1ᵧ = reshape(u1ᵧ, (length(mesh.centers[1])+1, length(mesh.centers[2])+1) )'
            reshaped_u2ₒ = reshape(u2ₒ, (length(mesh.centers[1])+1, length(mesh.centers[2])+1) )'
            reshaped_u2ᵧ = reshape(u2ᵧ, (length(mesh.centers[1])+1, length(mesh.centers[2])+1) )'
            
            ax1 = Axis(fig[1, 1], title="Diphasic Steady - Phase 1 - Bulk", xlabel="x", ylabel="y", aspect = DataAspect())
            hm1 = heatmap!(ax1, mesh.centers[1], mesh.centers[2], reshaped_u1ₒ, colormap=:viridis)
            cb1 = Colorbar(fig[1, 2], hm1, label="Phase 1 Bulk Temperature")
            #contour!(ax1, mesh.nodes[1], mesh.nodes[2], Z_sdf, levels=[0.0], color=:red, linewidth=2, label="SDF=0")

            ax2 = Axis(fig[1, 3], title="Diphasic Steady - Phase 1 - Interface", xlabel="x", ylabel="y", aspect = DataAspect())
            hm2 = heatmap!(ax2, mesh.centers[1], mesh.centers[2], reshaped_u1ᵧ, colormap=:viridis)
            cb2 = Colorbar(fig[1, 4], hm2, label="Phase 1 Interface Temperature")
            #contour!(ax2, mesh.nodes[1], mesh.nodes[2], Z_sdf, levels=[0.0], color=:red, linewidth=2, label="SDF=0")

            ax3 = Axis(fig[2, 1], title="Diphasic Steady - Phase 2 - Bulk", xlabel="x", ylabel="y", aspect = DataAspect())
            hm3 = heatmap!(ax3, mesh.centers[1], mesh.centers[2], reshaped_u2ₒ, colormap=:viridis)
            cb3 = Colorbar(fig[2, 2], hm3, label="Phase 2 Bulk Temperature")
            #contour!(ax3, mesh.nodes[1], mesh.nodes[2], Z_sdf, levels=[0.0], color=:red, linewidth=2, label="SDF=0")

            ax4 = Axis(fig[2, 3], title="Diphasic Steady - Phase 2 - Interface", xlabel="x", ylabel="y", aspect = DataAspect())
            hm4 = heatmap!(ax4, mesh.centers[1], mesh.centers[2], reshaped_u2ᵧ, colormap=:viridis)
            cb4 = Colorbar(fig[2, 4], hm4, label="Phase 2 Interface Temperature")
            #contour!(ax4, mesh.nodes[1], mesh.nodes[2], Z_sdf, levels=[0.0], color=:red, linewidth=2, label="SDF=0")

        end
        display(fig)
    else
        # Tracé unsteady en 2D
        if is_monophasic
            cell_types = capacity.cell_types
            uₒ = solver.states[state_i][1:length(solver.states[1]) ÷ 2]
            uᵧ = solver.states[state_i][length(solver.states[1]) ÷ 2 + 1:end]

            # Désactiver les cellules désactivées
            uₒ[cell_types .== 0] .= NaN
            uᵧ[cell_types .== 0] .= NaN
            uᵧ[cell_types .== 1] .= NaN

            reshaped_uₒ = reshape(uₒ, (length(mesh.centers[1])+1, length(mesh.centers[2])+1) )'
            reshaped_uᵧ = reshape(uᵧ, (length(mesh.centers[1])+1, length(mesh.centers[2])+1) )'

            fig = Figure(size=(800, 600))

            ax1 = Axis(fig[1, 1], title="Monophasic Unsteady Diffusion - Bulk", xlabel="x", ylabel="y", aspect = DataAspect())
            hm1 = heatmap!(ax1, mesh.centers[1], mesh.centers[2], reshaped_uₒ, colormap=:viridis)
            #contour!(ax1, mesh.nodes[1], mesh.nodes[2], Z_sdf, levels=[0.0], color=:red, linewidth=2, label="SDF=0")
            Colorbar(fig[1, 2], hm1, label="Bulk Temperature")

            ax2 = Axis(fig[1, 3], title="Monophasic Unsteady Diffusion - Interface", xlabel="x", ylabel="y", aspect = DataAspect())
            hm2 = heatmap!(ax2, mesh.centers[1], mesh.centers[2], reshaped_uᵧ, colormap=:viridis)
            #contour!(ax2, mesh.nodes[1], mesh.nodes[2], Z_sdf, levels=[0.0], color=:red, linewidth=2, label="SDF=0")
            Colorbar(fig[1, 4], hm2, label="Interface Temperature")

            display(fig)  
        else
            # Plot Last State
            cell_types = capacity.cell_types
            u1ₒ = solver.states[state_i][1:length(solver.states[1]) ÷ 4]
            u1ᵧ = solver.states[state_i][length(solver.states[1]) ÷ 4 + 1:2*length(solver.states[1]) ÷ 4]
            u2ₒ = solver.states[state_i][2*length(solver.states[1]) ÷ 4 + 1:3*length(solver.states[1]) ÷ 4]
            u2ᵧ = solver.states[state_i][3*length(solver.states[1]) ÷ 4 + 1:end]

            # Désactiver les cellules désactivées
            u1ₒ[cell_types .== 0] .= NaN
            u1ᵧ[cell_types .== 0] .= NaN
            u1ᵧ[cell_types .== 1] .= NaN
            u2ₒ[cell_types .== 1] .= NaN
            u2ᵧ[cell_types .== 1] .= NaN
            u2ᵧ[cell_types .== 0] .= NaN

            reshaped_u1ₒ = reshape(u1ₒ, (length(mesh.centers[1])+1, length(mesh.centers[2])+1) )'
            reshaped_u1ᵧ = reshape(u1ᵧ, (length(mesh.centers[1])+1, length(mesh.centers[2])+1) )'
            reshaped_u2ₒ = reshape(u2ₒ, (length(mesh.centers[1])+1, length(mesh.centers[2])+1) )'
            reshaped_u2ᵧ = reshape(u2ᵧ, (length(mesh.centers[1])+1, length(mesh.centers[2])+1) )'

            fig = Figure(size=(800, 600))

            x, y = range(mesh.x0[1], stop=mesh.x0[1]+mesh.h[1][1]*length(mesh.h[1]), length=length(mesh.h[1])+1), range(mesh.x0[2], stop=mesh.x0[2]+mesh.h[2][1]*length(mesh.h[2]), length=length(mesh.h[2])+1)

            ax1 = Axis(fig[1, 1], title="Diphasic Unsteady - Phase 1 - Bulk", xlabel="x", ylabel="y", aspect = DataAspect())
            hm1 = heatmap!(ax1, x, y, reshaped_u1ₒ, colormap=:viridis)
            Colorbar(fig[1, 2], hm1, label="Phase 1 Bulk Temperature")

            ax2 = Axis(fig[1, 3], title="Diphasic Unsteady - Phase 1 - Interface", xlabel="x", ylabel="y", aspect = DataAspect())
            hm2 = heatmap!(ax2, x, y, reshaped_u1ᵧ, colormap=:viridis)
            Colorbar(fig[1, 4], hm2, label="Phase 1 Interface Temperature")

            ax3 = Axis(fig[2, 1], title="Diphasic Unsteady - Phase 2 - Bulk", xlabel="x", ylabel="y", aspect = DataAspect())
            hm3 = heatmap!(ax3, x, y, reshaped_u2ₒ, colormap=:viridis)
            Colorbar(fig[2, 2], hm3, label="Phase 2 Bulk Temperature")

            ax4 = Axis(fig[2, 3], title="Diphasic Unsteady - Phase 2 - Interface", xlabel="x", ylabel="y", aspect = DataAspect())
            hm4 = heatmap!(ax4, x, y, reshaped_u2ᵧ, colormap=:viridis)
            Colorbar(fig[2, 4], hm4, label="Phase 2 Interface Temperature")

            display(fig)
        end
    end
end

function plot_solution(solver, mesh::CartesianMesh{3}, body::Body, capacity::Capacity; state_i=1)
    Z_sdf = [body.sdf(xi, yi, zi) for zi in mesh.nodes[3], yi in mesh.nodes[2], xi in mesh.nodes[1]]

    is_steady = solver.time_type == Steady
    is_monophasic = solver.phase_type == Monophasic

    # Tracé en 3D
    if is_steady
        if is_monophasic
            fig = Figure()
            ax = LScene(fig[1, 1], show_axis=false)

            nx, ny, nz = length(mesh.centers[1]), length(mesh.centers[2]), length(mesh.centers[3])
            x = LinRange(mesh.x0[1], mesh.x0[1]+mesh.h[1][1]*nx, nx+1)
            y = LinRange(mesh.x0[2], mesh.x0[2]+mesh.h[2][1]*ny, ny+1)
            z = LinRange(mesh.x0[3], mesh.x0[3]+mesh.h[3][1]*nz, nz+1)

            sgrid = SliderGrid(
                fig[2, 1],
                (label = "yz plane - x axis", range = 1:length(x)),
                (label = "xz plane - y axis", range = 1:length(y)),
                (label = "xy plane - z axis", range = 1:length(z)),
            )

            lo = sgrid.layout
            nc = ncols(lo)

            vol = solver.x[1:length(solver.x) ÷ 2]
            vol = reshape(vol, nx+1, ny+1, nz+1)
            plt = volumeslices!(ax, x, y, z, vol)
            Colorbar(fig[1, 2], plt)

            # connect sliders to `volumeslices` update methods
            sl_yz, sl_xz, sl_xy = sgrid.sliders

            on(sl_yz.value) do v; plt[:update_yz][](v) end
            on(sl_xz.value) do v; plt[:update_xz][](v) end
            on(sl_xy.value) do v; plt[:update_xy][](v) end

            set_close_to!(sl_yz, .5length(x))
            set_close_to!(sl_xz, .5length(y))
            set_close_to!(sl_xy, .5length(z))

            # add toggles to show/hide heatmaps
            hmaps = [plt[Symbol(:heatmap_, s)][] for s ∈ (:yz, :xz, :xy)]
            toggles = [Toggle(lo[i, nc + 1], active = true) for i ∈ 1:length(hmaps)]

            map(zip(hmaps, toggles)) do (h, t)
                connect!(h.visible, t.active)
            end

            # cam3d!(ax.scene, projectiontype=Makie.Orthographic)

            display(fig)
        else # Diphasic
            fig = Figure()
            ax1 = LScene(fig[1, 1], show_axis=false)
            ax2 = LScene(fig[1, 3], show_axis=false)

            nx, ny, nz = length(mesh.centers[1]), length(mesh.centers[2]), length(mesh.centers[3])
            x = LinRange(mesh.x0[1], mesh.x0[1]+mesh.h[1][1]*nx, nx+1)
            y = LinRange(mesh.x0[2], mesh.x0[2]+mesh.h[2][1]*ny, ny+1)
            z = LinRange(mesh.x0[3], mesh.x0[3]+mesh.h[3][1]*nz, nz+1)

            sgrid = SliderGrid(
                fig[2, 1],
                (label = "yz plane - x axis", range = 1:length(x)),
                (label = "xz plane - y axis", range = 1:length(y)),
                (label = "xy plane - z axis", range = 1:length(z)),
            )

            lo = sgrid.layout
            nc = ncols(lo)

            vol1 = solver.x[1:length(solver.x) ÷ 4]
            vol1 = reshape(vol1, nx+1, ny+1, nz+1)
            plt1 = volumeslices!(ax1, x, y, z, vol1)
            Colorbar(fig[1, 2], plt1)

            vol2 = solver.x[2*length(solver.x) ÷ 4 + 1:3*length(solver.x) ÷ 4]
            vol2 = reshape(vol2, nx+1, ny+1, nz+1)
            plt2 = volumeslices!(ax2, x, y, z, vol2)
            Colorbar(fig[1, 4], plt2)

            # connect sliders to `volumeslices` update methods
            sl_yz, sl_xz, sl_xy = sgrid.sliders

            on(sl_yz.value) do v
                plt1[:update_yz][](v) 
                plt2[:update_yz][](v)
            end
            on(sl_xz.value) do v
                plt1[:update_xz][](v)
                plt2[:update_xz][](v)
            end
            on(sl_xy.value) do v
                plt1[:update_xy][](v)
                plt2[:update_xy][](v)
            end

            set_close_to!(sl_yz, .5length(x))
            set_close_to!(sl_xz, .5length(y))
            set_close_to!(sl_xy, .5length(z))

            # add toggles to show/hide heatmaps
            hmaps = [plt1[Symbol(:heatmap_, s)][] for s ∈ (:yz, :xz, :xy)]
            toggles = [Toggle(lo[i, nc + 1], active = true) for i ∈ 1:length(hmaps)]

            map(zip(hmaps, toggles)) do (h, t)
                connect!(h.visible, t.active)
            end

            # cam3d!(ax1.scene, projectiontype=Makie.Orthographic)
            # cam3d!(ax2.scene, projectiontype=Makie.Orthographic)
            display(fig)
        end
    else # Unsteady
        if is_monophasic
            fig = Figure()
            ax = LScene(fig[1, 1], show_axis=false)

            nx, ny, nz = length(mesh.centers[1]), length(mesh.centers[2]), length(mesh.centers[3])
            x = LinRange(mesh.x0[1], mesh.x0[1]+mesh.h[1][1]*nx, nx+1)
            y = LinRange(mesh.x0[2], mesh.x0[2]+mesh.h[2][1]*ny, ny+1)
            z = LinRange(mesh.x0[3], mesh.x0[3]+mesh.h[3][1]*nz, nz+1)

            sgrid = SliderGrid(
                fig[2, 1],
                (label = "yz plane - x axis", range = 1:length(x)),
                (label = "xz plane - y axis", range = 1:length(y)),
                (label = "xy plane - z axis", range = 1:length(z)),
            )

            lo = sgrid.layout
            nc = ncols(lo)

            vol = solver.states[state_i][1:length(solver.states[state_i]) ÷ 2]
            vol = reshape(vol, nx+1, ny+1, nz+1)
            plt = volumeslices!(ax, x, y, z, vol)
            Colorbar(fig[1, 2], plt)

            # connect sliders to `volumeslices` update methods
            sl_yz, sl_xz, sl_xy = sgrid.sliders

            on(sl_yz.value) do v; plt[:update_yz][](v) end
            on(sl_xz.value) do v; plt[:update_xz][](v) end
            on(sl_xy.value) do v; plt[:update_xy][](v) end

            set_close_to!(sl_yz, .5length(x))
            set_close_to!(sl_xz, .5length(y))
            set_close_to!(sl_xy, .5length(z))

            # add toggles to show/hide heatmaps
            hmaps = [plt[Symbol(:heatmap_, s)][] for s ∈ (:yz, :xz, :xy)]
            toggles = [Toggle(lo[i, nc + 1], active = true) for i ∈ 1:length(hmaps)]

            map(zip(hmaps, toggles)) do (h, t)
                connect!(h.visible, t.active)
            end

            # cam3d!(ax.scene, projectiontype=Makie.Orthographic)  

            display(fig)

        else # Diphasic
            fig = Figure()
            ax1 = LScene(fig[1, 1], show_axis=false)
            ax2 = LScene(fig[1, 3], show_axis=false)

            nx, ny, nz = length(mesh.centers[1]), length(mesh.centers[2]), length(mesh.centers[3])
            x = LinRange(mesh.x0[1], mesh.x0[1]+mesh.h[1][1]*nx, nx+1)
            y = LinRange(mesh.x0[2], mesh.x0[2]+mesh.h[2][1]*ny, ny+1)
            z = LinRange(mesh.x0[3], mesh.x0[3]+mesh.h[3][1]*nz, nz+1)

            sgrid = SliderGrid(
                fig[2, 1],
                (label = "yz plane - x axis", range = 1:length(x)),
                (label = "xz plane - y axis", range = 1:length(y)),
                (label = "xy plane - z axis", range = 1:length(z)),
            )

            lo = sgrid.layout
            nc = ncols(lo)

            vol1 = solver.states[state_i][1:length(solver.states[state_i]) ÷ 4]
            vol1 = reshape(vol1, nx+1, ny+1, nz+1)
            plt1 = volumeslices!(ax1, x, y, z, vol1)
            Colorbar(fig[1, 2], plt1)

            vol2 = solver.states[state_i][2*length(solver.states[state_i]) ÷ 4 + 1:3*length(solver.states[state_i]) ÷ 4]
            vol2 = reshape(vol2, nx+1, ny+1, nz+1)
            plt2 = volumeslices!(ax2, x, y, z, vol2)
            Colorbar(fig[1, 4], plt2)

            # connect sliders to `volumeslices` update methods
            sl_yz, sl_xz, sl_xy = sgrid.sliders

            on(sl_yz.value) do v
                plt1[:update_yz][](v) 
                plt2[:update_yz][](v)
            end
            on(sl_xz.value) do v
                plt1[:update_xz][](v)
                plt2[:update_xz][](v)
            end
            on(sl_xy.value) do v
                plt1[:update_xy][](v)
                plt2[:update_xy][](v)
            end

            set_close_to!(sl_yz, .5length(x))
            set_close_to!(sl_xz, .5length(y))
            set_close_to!(sl_xy, .5length(z))

            # add toggles to show/hide heatmaps
            hmaps = [plt1[Symbol(:heatmap_, s)][] for s ∈ (:yz, :xz, :xy)]
            toggles = [Toggle(lo[i, nc + 1], active = true) for i ∈ 1:length(hmaps)]

            map(zip(hmaps, toggles)) do (h, t)
                connect!(h.visible, t.active)
            end

            # cam3d!(ax1.scene, projectiontype=Makie.Orthographic)
            # cam3d!(ax2.scene, projectiontype=Makie.Orthographic)
            display(fig)
        end
    end
end


function plot_profile(solver, nx, ny, x0, lx, y0, ly, x; ntime=4)
    # Récupérer les coordonnées y
    y = range(y0, stop=ly, length=ny+1)

    # Récupérer les indices x les plus proches
    x_idx = round(Int, (x - x0) / (lx / nx))

    # Sélectionner 4 instants équidistants
    total_states = length(solver.states)
    indices = round.(Int, LinRange(5, total_states, ntime))

    fig = Figure(size = (800, 600))

    # Créer les sous-graphiques pour les deux phases
    ax1 = Axis(fig[1, 1], title = "Phase 1 - Bulk", xlabel = "y", ylabel = "u1ₒ")
    ax2 = Axis(fig[2, 1], title = "Phase 2 - Bulk", xlabel = "y", ylabel = "u2ₒ")

    for (i, idx) in enumerate(indices)
        # Récupérer les valeurs de la solution pour les 4 phases à l'instant idx
        state = solver.states[idx]
        u1ₒ = reshape(state[1:length(state) ÷ 4], (nx + 1, ny + 1))[x_idx, :]
        u2ₒ = reshape(state[2*length(state) ÷ 4 + 1:3*length(state) ÷ 4], (nx + 1, ny + 1))[x_idx, :]

        # Tracer les lignes pour chaque instant
        lines!(ax1, y, u1ₒ, label = "Instant $idx")
        lines!(ax2, y, u2ₒ, label = "Instant $idx")
    end

    # Ajouter des légendes
    axislegend(ax1, position = :rt)
    axislegend(ax2, position = :rt)

    # Afficher le graphique
    display(fig)
end

function plot_profile(solver::Solver, mesh::CartesianMesh; x=1.0)
    # Plot Profile for 2D or 3D

    dims = length(mesh.centers)

    if dims == 2
        # Récupérer les coordonnées y
        y = range(mesh.x0[2], stop=mesh.x0[2]+mesh.h[2][1]*length(mesh.h[2]), length=length(mesh.h[2])+1)

        # Récupérer les indices x les plus proches
        x_idx = round(Int, (x - mesh.x0[1]) / mesh.h[1][1])

        # Sélectionner 4 instants équidistants
        total_states = length(solver.states)
        indices = round.(Int, LinRange(5, total_states, 4))

        fig = Figure(size = (800, 600))

        # Créer les sous-graphiques pour les deux phases
        ax1 = Axis(fig[1, 1], title = "Phase 1 - Bulk", xlabel = "y", ylabel = "u1ₒ")
        ax2 = Axis(fig[2, 1], title = "Phase 2 - Bulk", xlabel = "y", ylabel = "u2ₒ")

        for (i, idx) in enumerate(indices)
            # Récupérer les valeurs de la solution pour les 4 phases à l'instant idx
            state = solver.states[idx]
            u1ₒ = reshape(state[1:length(state) ÷ 4], (length(mesh.centers[1])+1, length(mesh.centers[2])+1))[x_idx, 1:div(length(mesh.centers[2])+1, 2)]
            u2ₒ = reshape(state[2*length(state) ÷ 4 + 1:3*length(state) ÷ 4], (length(mesh.centers[1])+1, length(mesh.centers[2])+1))[x_idx, 1:div(length(mesh.centers[2])+1, 2)]

            # Tracer les lignes pour chaque instant
            lines!(ax1, y[1:div(length(y), 2)], u1ₒ, label = "Instant $idx")
            lines!(ax2, y[1:div(length(y), 2)], u2ₒ, label = "Instant $idx")
            # Ajouter des points pour u10 et u20 (10% des points)
            scatter!(ax1, y[1:5:div(length(y), 2)], u1ₒ[1:5:div(length(y), 2)], color=:red, marker=:circle, label="u10 Instant $idx")
            scatter!(ax2, y[1:5:div(length(y), 2)], u2ₒ[1:5:div(length(y), 2)], color=:blue, marker=:circle, label="u20 Instant $idx")
        end

        # Ajouter des légendes
        axislegend(ax1, position = :rt)
        axislegend(ax2, position = :rt)

        # Afficher le graphique
        display(fig)

    elseif dims == 3
        println("Plotting Profile for 3D is not supported yet.")

    else
        error("Dimension non supportée: $dims")
    end
end


# Plot the mesh
function plot_mesh(mesh_p, mesh_u)
    fig = Figure(resolution = (800, 800))
    ax = Axis(fig[1, 1], title = "Mesh Plot", xlabel = "x", ylabel = "y", aspect = DataAspect())

    # Plot mesh_p centers
    scatter!(ax, repeat(mesh_p.centers[1], inner=length(mesh_p.centers[2])), repeat(mesh_p.centers[2], outer=length(mesh_p.centers[1])), color = :blue, label = "mesh_p", markersize=15, marker = 'o')

    # Plot mesh_u centers triangles
    scatter!(ax, repeat(mesh_u.centers[1], inner=length(mesh_u.centers[2])), repeat(mesh_u.centers[2], outer=length(mesh_u.centers[1])), color = :red, label = "mesh_u", markersize=15, marker = :rtriangle)

    axislegend(ax)
    display(fig)
end


function plot_mesh(mesh_p, mesh_u, mesh_v)
    fig = Figure(resolution = (800, 800))
    ax = Axis(fig[1, 1], title = "Mesh Plot", xlabel = "x", ylabel = "y", aspect = DataAspect())

    # Plot mesh_p centers
    scatter!(ax, repeat(mesh_p.centers[1], inner=length(mesh_p.centers[2])), repeat(mesh_p.centers[2], outer=length(mesh_p.centers[1])), color = :blue, label = "mesh_p", markersize=15, marker = 'o')

    # Plot mesh_u centers triangles
    scatter!(ax, repeat(mesh_u.centers[1], inner=length(mesh_u.centers[2])), repeat(mesh_u.centers[2], outer=length(mesh_u.centers[1])), color = :red, label = "mesh_u", markersize=15, marker = :rtriangle)

    # Plot mesh_v centers
    scatter!(ax, repeat(mesh_v.centers[1], inner=length(mesh_v.centers[2])), repeat(mesh_v.centers[2], outer=length(mesh_v.centers[1])), color = :green, label = "mesh_v", markersize=15, marker = :utriangle)

    axislegend(ax)
    display(fig)
end

function plot_mesh(mesh_p, mesh_u, mesh_v, mesh_w)
    fig = Figure(resolution = (800, 800))
    ax = Axis3(fig[1, 1], title = "Mesh Plot", xlabel = "x", ylabel = "y", zlabel = "z", aspect = DataAspect())

    # Plot mesh_p centers
    scatter!(ax, repeat(mesh_p.centers[1], inner=length(mesh_p.centers[2])), repeat(mesh_p.centers[2], outer=length(mesh_p.centers[1])), repeat(mesh_p.centers[3], outer=length(mesh_p.centers[1])*length(mesh_p.centers[2])), color = :blue, label = "mesh_p", markersize=15, marker = 'o')

    # Plot mesh_u centers triangles
    scatter!(ax, repeat(mesh_u.centers[1], inner=length(mesh_u.centers[2])), repeat(mesh_u.centers[2], outer=length(mesh_u.centers[1])), repeat(mesh_u.centers[3], outer=length(mesh_u.centers[1])*length(mesh_u.centers[2])), color = :red, label = "mesh_u", markersize=15, marker = :rtriangle)

    # Plot mesh_v centers
    scatter!(ax, repeat(mesh_v.centers[1], inner=length(mesh_v.centers[2])), repeat(mesh_v.centers[2], outer=length(mesh_v.centers[1])), repeat(mesh_v.centers[3], outer=length(mesh_v.centers[1])*length(mesh_v.centers[2])), color = :green, label = "mesh_v", markersize=15, marker = :utriangle)

    # Plot mesh_w centers
    scatter!(ax, repeat(mesh_w.centers[1], inner=length(mesh_w.centers[2])), repeat(mesh_w.centers[2], outer=length(mesh_w.centers[1])), repeat(mesh_w.centers[3], outer=length(mesh_w.centers[1])*length(mesh_w.centers[2])), color = :yellow, label = "mesh_w", markersize=15, marker = :dtriangle)

    axislegend(ax)
    display(fig)
end