"""
    write_vtk(filename::String, mesh::CartesianMesh, solver::Solver)

Writes a VTK file based on the solver's time, phase, and equation types.

# Arguments
- `filename::String` : Name of the VTK file to create.
- `mesh::CartesianMesh` : Mesh used in the simulation.
- `solver::Solver` : Simulation results containing the data to write.
"""
function write_vtk(filename::String, mesh::CartesianMesh, solver::Solver)
    if length(mesh.centers) == 1
        vtk_grid(filename, mesh.centers[1], ones(length(mesh.centers[1]))) do vtk
            if solver.time_type == Steady && solver.phase_type == Monophasic && solver.equation_type == Diffusion
                # Cas Steady, Monophasic, Diffusion
                vtk["Temperature_b"] = solver.x[1:div(length(solver.x), 2)]
                vtk["Temperature_g"] = solver.x[div(length(solver.x), 2) + 1:end]
                println("Fichier VTK steady monophasic écrit : $filename.vtr")

            elseif solver.time_type == Steady && solver.phase_type == Diphasic && solver.equation_type == Diffusion
                # Cas Steady, Diphasic, Diffusion
                part = div(length(solver.x), 4)
                vtk["Temperature_1_b"] = solver.x[1:part]
                vtk["Temperature_1_g"] = solver.x[part + 1:2 * part]
                vtk["Temperature_2_b"] = solver.x[2 * part + 1:3 * part]
                vtk["Temperature_2_g"] = solver.x[3 * part + 1:end]
                println("Fichier VTK steady diphasic écrit : $filename.vtr")

            elseif solver.time_type == Unsteady && solver.phase_type == Monophasic && solver.equation_type == Diffusion
                # Cas Unsteady, Monophasic, Diffusion
                for (i, state) in enumerate(solver.states)
                    vtk["Temperature_b_$i"] = state[1:div(length(state), 2)]
                    vtk["Temperature_g_$i"] = state[div(length(state), 2) + 1:end]
                end
                println("Fichier VTK unsteady monophasic écrit : $filename.vtr")

            elseif solver.time_type == Unsteady && solver.phase_type == Diphasic && solver.equation_type == Diffusion
                # Cas Unsteady, Diphasic, Diffusion
                part = div(length(solver.x), 4)
                for (i, state) in enumerate(solver.states)
                    vtk["Temperature_1_b_$i"] = state[1:part]
                    vtk["Temperature_1_g_$i"] = state[part + 1:2 * part]
                    vtk["Temperature_2_b_$i"] = state[2 * part + 1:3 * part]
                    vtk["Temperature_2_g_$i"] = state[3 * part + 1:end]
                end
                println("Fichier VTK unsteady diphasic écrit : $filename.vtr")

            # Ajoutez d'autres combinaisons si nécessaire
            else
                error("La combinaison de TimeType, PhaseType et EquationType n'est pas supportée.")
            end
        end
    elseif length(mesh.centers) == 2
        vtk_grid(filename, mesh.centers[1], mesh.centers[2]) do vtk
            if solver.time_type == Steady && solver.phase_type == Monophasic && solver.equation_type == Diffusion
                # Cas Steady, Monophasic, Diffusion
                vtk["Temperature_b"] = solver.x[1:div(length(solver.x), 2)]
                vtk["Temperature_g"] = solver.x[div(length(solver.x), 2) + 1:end]
                println("Fichier VTK steady monophasic écrit : $filename.vtr")
            
            elseif solver.time_type == Steady && solver.phase_type == Diphasic && solver.equation_type == Diffusion
                # Cas Steady, Diphasic, Diffusion
                part = div(length(solver.x), 4)
                vtk["Temperature_1_b"] = solver.x[1:part]
                vtk["Temperature_1_g"] = solver.x[part + 1:2 * part]
                vtk["Temperature_2_b"] = solver.x[2 * part + 1:3 * part]
                vtk["Temperature_2_g"] = solver.x[3 * part + 1:end]
                println("Fichier VTK steady diphasic écrit : $filename.vtr")
            
            elseif solver.time_type == Unsteady && solver.phase_type == Monophasic && solver.equation_type == Diffusion
                # Cas Unsteady, Monophasic, Diffusion
                for (i, state) in enumerate(solver.states)
                    vtk["Temperature_b_$i"] = state[1:div(length(state), 2)]
                    vtk["Temperature_g_$i"] = state[div(length(state), 2) + 1:end]
                end
                println("Fichier VTK unsteady monophasic écrit : $filename.vtr")
            
            elseif solver.time_type == Unsteady && solver.phase_type == Diphasic && solver.equation_type == Diffusion
                # Cas Unsteady, Diphasic, Diffusion
                part = div(length(solver.x), 4)
                for (i, state) in enumerate(solver.states)
                    vtk["Temperature_1_b_$i"] = state[1:part]
                    vtk["Temperature_1_g_$i"] = state[part + 1:2 * part]
                    vtk["Temperature_2_b_$i"] = state[2 * part + 1:3 * part]
                    vtk["Temperature_2_g_$i"] = state[3 * part + 1:end]
                end
                println("Fichier VTK unsteady diphasic écrit : $filename.vtr")
    
            # Ajoutez d'autres combinaisons si nécessaire
            else
                error("La combinaison de TimeType, PhaseType et EquationType n'est pas supportée.")
            end
        end
    elseif length(mesh.centers) == 3
        vtk_grid(filename, mesh.centers[1], mesh.centers[2], mesh.centers[3]) do vtk
            # Code for 3D case
            # ...
        end
    else
        error("Invalid number of dimensions for mesh.centers.")
    end
end
