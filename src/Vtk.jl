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
        if solver.time_type == Steady && solver.phase_type == Monophasic
            # Steady, Monophasic, Diffusion
            vtk_grid(filename, 0:1:length(mesh.centers[1])) do vtk
                vtk["Temperature_b"] = solver.x[1:length(solver.x) ÷ 2]
                vtk["Temperature_g"] = solver.x[length(solver.x) ÷ 2 + 1:end]
                println("VTK file written : $filename.vti")
            end
        elseif solver.time_type == Steady && solver.phase_type == Diphasic
            # Steady, Diphasic, Diffusion
            part = div(length(solver.x), 4)
            vtk_grid(filename, 0:1:length(mesh.centers[1])) do vtk
                vtk["Temperature_1_b"] = solver.x[1:part]
                vtk["Temperature_1_g"] = solver.x[part + 1:2 * part]
                vtk["Temperature_2_b"] = solver.x[2 * part + 1:3 * part]
                vtk["Temperature_2_g"] = solver.x[3 * part + 1:end]
                println("VTK file written : $filename.vti")
            end
        elseif solver.time_type == Unsteady && solver.phase_type == Monophasic
            # Unsteady, Monophasic, Diffusion
            pvd = paraview_collection(filename)
            for (i, state) in enumerate(solver.states)
                vtk_grid(filename * "_$i", 0:1:length(mesh.centers[1])) do vtk
                    vtk["Temperature_b"] = state[1:length(state) ÷ 2]
                    vtk["Temperature_g"] = state[length(state) ÷ 2 + 1:end]
                    pvd[i] = vtk
                end
            end
            vtk_save(pvd)
            println("VTK file written : $filename.pvd")
        elseif solver.time_type == Unsteady && solver.phase_type == Diphasic
            # Unsteady, Diphasic, Diffusion
            pvd = paraview_collection(filename)
            part = div(length(solver.x), 4)
            for (i, state) in enumerate(solver.states)
                vtk_grid(filename * "_$i", 0:1:length(mesh.centers[1])) do vtk
                    vtk["Temperature_1_b"] = state[1:part]
                    vtk["Temperature_1_g"] = state[part + 1:2 * part]
                    vtk["Temperature_2_b"] = state[2 * part + 1:3 * part]
                    vtk["Temperature_2_g"] = state[3 * part + 1:end]
                    pvd[i] = vtk
                end
            end
            vtk_save(pvd)
            println("VTK file written : $filename.pvd")
        else
            error("Combination of TimeType, PhaseType, and EquationType not supported.")
        end
    elseif length(mesh.centers) == 2
        if solver.time_type == Steady && solver.phase_type == Monophasic 
            # Cas Steady, Monophasic, Diffusion
            vtk_grid(filename, 0:1:length(mesh.centers[1]), 0:1:length(mesh.centers[2])) do vtk
                vtk["Temperature_b"] = reshape(solver.x[1:length(solver.x) ÷ 2], length(mesh.centers[1])+1, length(mesh.centers[2])+1)
                vtk["Temperature_g"] = reshape(solver.x[length(solver.x) ÷ 2 + 1:end], length(mesh.centers[1])+1, length(mesh.centers[2])+1)
                println("Fichier VTK steady monophasic écrit : $filename.vti")
            end
        elseif solver.time_type == Steady && solver.phase_type == Diphasic 
            # Cas Steady, Diphasic, Diffusion
            part = div(length(solver.x), 4)
            vtk_grid(filename, 0:1:length(mesh.centers[1]), 0:1:length(mesh.centers[2])) do vtk
                vtk["Temperature_1_b"] = reshape(solver.x[1:part], length(mesh.centers[1])+1, length(mesh.centers[2])+1)
                vtk["Temperature_1_g"] = reshape(solver.x[part + 1:2 * part], length(mesh.centers[1])+1, length(mesh.centers[2])+1)
                vtk["Temperature_2_b"] = reshape(solver.x[2 * part + 1:3 * part], length(mesh.centers[1])+1, length(mesh.centers[2])+1)
                vtk["Temperature_2_g"] = reshape(solver.x[3 * part + 1:end], length(mesh.centers[1])+1, length(mesh.centers[2])+1)
                println("Fichier VTK steady diphasic écrit : $filename.vti")
            end
        elseif solver.time_type == Unsteady && solver.phase_type == Monophasic 
            pvd = paraview_collection(filename)
            # Cas Unsteady, Monophasic, Diffusion
            for (i, state) in enumerate(solver.states)
                vtk_grid(filename * "_$i", 0:1:length(mesh.centers[1]), 0:1:length(mesh.centers[2])) do vtk
                    vtk["Temperature_b"] = reshape(state[1:length(state) ÷ 2], length(mesh.centers[1])+1, length(mesh.centers[2])+1)
                    vtk["Temperature_g"] = reshape(state[length(state) ÷ 2 + 1:end], length(mesh.centers[1])+1, length(mesh.centers[2])+1)
                    pvd[i] = vtk
                end
            end
            vtk_save(pvd)
            println("Fichier VTK unsteady monophasic écrit : $filename.pvd")
        elseif solver.time_type == Unsteady && solver.phase_type == Diphasic 
            pvd = paraview_collection(filename)
            # Cas Unsteady, Diphasic, Diffusion
            part = div(length(solver.x), 4)
            for (i, state) in enumerate(solver.states)
                vtk_grid(filename * "_$i", 0:1:length(mesh.centers[1]), 0:1:length(mesh.centers[2])) do vtk
                    vtk["Temperature_1_b"] = reshape(state[1:part], length(mesh.centers[1])+1, length(mesh.centers[2])+1)
                    vtk["Temperature_1_g"] = reshape(state[part + 1:2 * part], length(mesh.centers[1])+1, length(mesh.centers[2])+1)
                    vtk["Temperature_2_b"] = reshape(state[2 * part + 1:3 * part], length(mesh.centers[1])+1, length(mesh.centers[2])+1)
                    vtk["Temperature_2_g"] = reshape(state[3 * part + 1:end], length(mesh.centers[1])+1, length(mesh.centers[2])+1)
                    pvd[i] = vtk
                end
            end
            vtk_save(pvd)
            println("Fichier VTK unsteady diphasic écrit : $filename.pvd")            
        else
            error("La combinaison de TimeType, PhaseType et EquationType n'est pas supportée.")
        end
    elseif length(mesh.centers) == 3
        if solver.time_type == Steady && solver.phase_type == Monophasic 
            # Cas Steady, Monophasic, Diffusion
            vtk_grid(filename, 0:1:length(mesh.centers[1]), 0:1:length(mesh.centers[2]), 0:1:length(mesh.centers[3])) do vtk
                vtk["Temperature_b"] = reshape(solver.x[1:length(solver.x) ÷ 2], length(mesh.centers[1])+1, length(mesh.centers[2])+1, length(mesh.centers[3])+1)
                vtk["Temperature_g"] = reshape(solver.x[length(solver.x) ÷ 2 + 1:end], length(mesh.centers[1])+1, length(mesh.centers[2])+1, length(mesh.centers[3])+1)
                println("Fichier VTK steady monophasic écrit : $filename.vti")
            end
        elseif solver.time_type == Steady && solver.phase_type == Diphasic 
            # Cas Steady, Diphasic, Diffusion
            part = div(length(solver.x), 4)
            vtk_grid(filename, 0:1:length(mesh.centers[1]), 0:1:length(mesh.centers[2]), 0:1:length(mesh.centers[3])) do vtk
                vtk["Temperature_1_b"] = reshape(solver.x[1:part], length(mesh.centers[1])+1, length(mesh.centers[2])+1, length(mesh.centers[3])+1)
                vtk["Temperature_1_g"] = reshape(solver.x[part + 1:2 * part], length(mesh.centers[1])+1, length(mesh.centers[2])+1, length(mesh.centers[3])+1)
                vtk["Temperature_2_b"] = reshape(solver.x[2 * part + 1:3 * part], length(mesh.centers[1])+1, length(mesh.centers[2])+1, length(mesh.centers[3])+1)
                vtk["Temperature_2_g"] = reshape(solver.x[3 * part + 1:end], length(mesh.centers[1])+1, length(mesh.centers[2])+1, length(mesh.centers[3])+1)
                println("Fichier VTK steady diphasic écrit : $filename.vti")
            end
        elseif solver.time_type == Unsteady && solver.phase_type == Monophasic 
            pvd = paraview_collection(filename)
            # Cas Unsteady, Monophasic, Diffusion
            for (i, state) in enumerate(solver.states)
                vtk_grid(filename * "_$i", 0:1:length(mesh.centers[1]), 0:1:length(mesh.centers[2]), 0:1:length(mesh.centers[3])) do vtk
                    vtk["Temperature_b"] = reshape(state[1:length(state) ÷ 2], length(mesh.centers[1])+1, length(mesh.centers[2])+1, length(mesh.centers[3])+1)
                    vtk["Temperature_g"] = reshape(state[length(state) ÷ 2 + 1:end], length(mesh.centers[1])+1, length(mesh.centers[2])+1, length(mesh.centers[3])+1)
                    pvd[i] = vtk
                end
            end
            vtk_save(pvd)
            println("Fichier VTK unsteady monophasic écrit : $filename.pvd")
        elseif solver.time_type == Unsteady && solver.phase_type == Diphasic 
            pvd = paraview_collection(filename)
            # Cas Unsteady, Diphasic, Diffusion
            part = div(length(solver.x), 4)
            for (i, state) in enumerate(solver.states)
                vtk_grid(filename * "_$i", 0:1:length(mesh.centers[1]), 0:1:length(mesh.centers[2]), 0:1:length(mesh.centers[3])) do vtk
                    vtk["Temperature_1_b"] = reshape(state[1:part], length(mesh.centers[1])+1, length(mesh.centers[2])+1, length(mesh.centers[3])+1)
                    vtk["Temperature_1_g"] = reshape(state[part + 1:2 * part], length(mesh.centers[1])+1, length(mesh.centers[2])+1, length(mesh.centers[3])+1)
                    vtk["Temperature_2_b"] = reshape(state[2 * part + 1:3 * part], length(mesh.centers[1])+1, length(mesh.centers[2])+1, length(mesh.centers[3])+1)
                    vtk["Temperature_2_g"] = reshape(state[3 * part + 1:end], length(mesh.centers[1])+1, length(mesh.centers[2])+1, length(mesh.centers[3])+1)
                    pvd[i] = vtk
                end
            end
            vtk_save(pvd)
            println("Fichier VTK unsteady diphasic écrit : $filename.pvd")
        else 
            error("La combinaison de TimeType, PhaseType et EquationType n'est pas supportée.")
        end
    else
        error("Invalid number of dimensions for mesh.centers.")
    end
end

function write_vtk(filename::String, mesh::CartesianMesh, solver::SolverVec)
    if length(mesh.centers) == 1
        if solver.time_type == Unsteady && solver.phase_type == Monophasic
            # Unsteady, Monophasic, Diffusion
            pvd = paraview_collection(filename)
            for (i, state) in enumerate(solver.states)
                vtk_grid(filename * "_$i", 0:1:length(mesh.centers[1])) do vtk
                    vtk["Velocity_x_o"] = state[1:length(state) ÷ 2]
                    vtk["Velocity_x_γ"] = state[length(state) ÷ 2 + 1:end]
                    pvd[i] = vtk
                end
            end
            vtk_save(pvd)
            println("VTK file written : $filename.pvd")
        else
            error("Combination of TimeType, PhaseType, and EquationType not supported.")
        end
    elseif length(mesh.centers) == 2
        if solver.time_type == Unsteady && solver.phase_type == Monophasic
            pvd_u = paraview_collection(filename * "_u")
            pvd_v = paraview_collection(filename * "_v")
            # Cas Unsteady, Monophasic, Diffusion
            for (i, state) in enumerate(solver.states)
                len_ux = length(mesh.centers[1]) * (length(mesh.centers[2]) + 1)
                len_uy = (length(mesh.centers[1]) + 1) * length(mesh.centers[2])

                uxₒ = state[1:len_ux]
                uxᵧ = state[len_ux + 1:2 * len_ux]
                uyₒ = state[2 * len_ux + 1:2 * len_ux + len_uy]
                uyᵧ = state[2 * len_ux + len_uy + 1:end]

                Uxₒ = reshape(uxₒ, length(mesh.centers[1]), length(mesh.centers[2]) + 1)
                Uxᵧ = reshape(uxᵧ, length(mesh.centers[1]), length(mesh.centers[2]) + 1)
                Uyₒ = reshape(uyₒ, length(mesh.centers[1]) + 1, length(mesh.centers[2]))
                Uyᵧ = reshape(uyᵧ, length(mesh.centers[1]) + 1, length(mesh.centers[2]))
                
                vtk_grid(filename * "_u_$i", 1:1:length(mesh.centers[1]), 0:1:length(mesh.centers[2])) do vtk_u
                    vtk_u["Velocity_x_o"] = Uxₒ
                    vtk_u["Velocity_x_γ"] = Uxᵧ
                    pvd_u[i] = vtk_u
                end
                
                vtk_grid(filename * "_v_$i", 0:1:length(mesh.centers[1]), 1:1:length(mesh.centers[2])) do vtk_v
                    vtk_v["Velocity_y_o"] = Uyₒ
                    vtk_v["Velocity_y_γ"] = Uyᵧ
                    pvd_v[i] = vtk_v
                end
            end
            vtk_save(pvd_u)
            vtk_save(pvd_v)
            println("VTK files written : filename_u.pvd and filename_v.pvd")
        else
            error("Combination of TimeType, PhaseType, and EquationType not supported.")
        end
    elseif length(mesh.centers) == 3
        if solver.time_type == Unsteady && solver.phase_type == Monophasic
            pvd_u = paraview_collection(filename * "_u")
            pvd_v = paraview_collection(filename * "_v")
            pvd_w = paraview_collection(filename * "_w")
            # Cas Unsteady, Monophasic, Diffusion
            for (i, state) in enumerate(solver.states)
                len_ux = length(mesh.centers[1]) * (length(mesh.centers[2]) + 1)
                len_uy = (length(mesh.centers[1]) + 1) * length(mesh.centers[2])
                len_uz = (length(mesh.centers[1]) + 1) * (length(mesh.centers[2]) + 1) * length(mesh.centers[3])

                uxₒ = state[1:len_ux]
                uxᵧ = state[len_ux + 1:2 * len_ux]
                uyₒ = state[2 * len_ux + 1:2 * len_ux + len_uy]
                uyᵧ = state[2 * len_ux + len_uy + 1:2 * len_ux + len_uy + len_uz]
                uzₒ = state[2 * len_ux + len_uy + len_uz + 1:2 * len_ux + len_uy + 2 * len_uz]
                uzᵧ = state[2 * len_ux + len_uy + 2 * len_uz + 1:end]


                Uxₒ = reshape(uxₒ, length(mesh.centers[1]), length(mesh.centers[2]) + 1, length(mesh.centers[3]))
                Uxᵧ = reshape(uxᵧ, length(mesh.centers[1]), length(mesh.centers[2]) + 1, length(mesh.centers[3]))
                Uyₒ = reshape(uyₒ, length(mesh.centers[1]) + 1, length(mesh.centers[2]), length(mesh.centers[3]))
                Uyᵧ = reshape(uyᵧ, length(mesh.centers[1]) + 1, length(mesh.centers[2]), length(mesh.centers[3]))
                Uzₒ = reshape(uzₒ, length(mesh.centers[1]) + 1, length(mesh.centers[2]) + 1, length(mesh.centers[3]))
                Uzᵧ = reshape(uzᵧ, length(mesh.centers[1]) + 1, length(mesh.centers[2]) + 1, length(mesh.centers[3]))

                vtk_u = vtk_grid(filename * "_u_$i", 1:1:length(mesh.centers[1]), 0:1:length(mesh.centers[2]), 0:1:length(mesh.centers[3])) do vtk
                    vtk["Velocity_x_o"] = Uxₒ
                    vtk["Velocity_x_γ"] = Uxᵧ
                    pvd_u[i] = vtk
                end

                vtk_v = vtk_grid(filename * "_v_$i", 0:1:length(mesh.centers[1]), 1:1:length(mesh.centers[2]), 0:1:length(mesh.centers[3])) do vtk
                    vtk["Velocity_y_o"] = Uyₒ
                    vtk["Velocity_y_γ"] = Uyᵧ
                    pvd_v[i] = vtk
                end

                vtk_w = vtk_grid(filename * "_w_$i", 0:1:length(mesh.centers[1]), 0:1:length(mesh.centers[2]), 1:1:length(mesh.centers[3])) do vtk
                    vtk["Velocity_z_o"] = Uzₒ
                    vtk["Velocity_z_γ"] = Uzᵧ
                    pvd_w[i] = vtk
                end
            end
            vtk_save(pvd_u)
            vtk_save(pvd_v)
            vtk_save(pvd_w)
            println("VTK files written : filename_u.pvd, filename_v.pvd, and filename_w.pvd")
        else
            error("Combination of TimeType, PhaseType, and EquationType not supported.")
        end
    else
        error("Invalid number of dimensions for mesh.centers.")
    end
end
