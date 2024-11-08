using Fliwer
using IterativeSolvers
using CairoMakie

# Fonction pour configurer et résoudre le problème pour un maillage donné
function solve_diffusion_unsteady(nx, ny, lx, ly, x0, y0, radius, center, Δt, Tend, u0)

    # Définir le domaine et le maillage
    domain = ((x0, lx), (y0, ly))
    mesh = CartesianMesh((nx, ny), (lx, ly), (x0, y0))

    # Définir le corps (disque)
    circle = Body((x, y, _=0) -> sqrt((x - center[1])^2 + (y - center[2])^2) - radius,
                  (x, y, _) -> (x, y), domain, false)
    circle_c = Body((x, y, _=0) -> -(sqrt((x - center[1])^2 + (y - center[2])^2) - radius),
                    (x, y, _) -> (x, y), domain, false)

    # Identifier les cellules contenant le corps
    identify!(mesh, circle)

    # Définir les capacités
    capacity = Capacity(circle, mesh)
    capacity_c = Capacity(circle_c, mesh)

    # Définir les opérateurs de diffusion
    operator = DiffusionOps(capacity.A, capacity.B, capacity.V, capacity.W, (nx + 1, ny + 1))
    operator_c = DiffusionOps(capacity_c.A, capacity_c.B, capacity_c.V, capacity_c.W, (nx + 1, ny + 1))

    # Définir les conditions aux limites
    bc = Dirichlet(0.0)
    bc_b = BorderConditions(Dict{Symbol, AbstractBoundary}(:left => bc, :right => bc, :top => bc, :bottom => bc))

    # Définir les conditions d'interface
    ic = InterfaceConditions(ScalarJump(1.0, 2.0, 0.0), FluxJump(1.0, 1.0, 0.0))

    # Définir le terme source
    f1 = (x, y, z, t) -> 0.0
    f2 = (x, y, z, t) -> 0.0

    # Définir les phases
    Fluide_1 = Phase(capacity, operator, f1, 1.0)
    Fluide_2 = Phase(capacity_c, operator_c, f2, 1.0)

    # Initial condition
    u0ₒ1 = ones((nx + 1) * (ny + 1))
    u0ᵧ1 = zeros((nx + 1) * (ny + 1))
    u0ₒ2 = zeros((nx + 1) * (ny + 1))
    u0ᵧ2 = zeros((nx + 1) * (ny + 1))
    u0_full = vcat(u0ₒ1, u0ᵧ1, u0ₒ2, u0ᵧ2)

    # Définir le solver
    solver = DiffusionUnsteadyDiph(Fluide_1, Fluide_2, bc_b, ic, Δt, Tend, u0_full)

    # Résoudre le problème
    solve!(solver, Fluide_1, Fluide_2, u0_full, Δt, Tend, bc_b, ic;
           method=IterativeSolvers.gmres,
           restart=100,
           maxiter=1000,
           verbose=false)

    return mesh, solver
end

# Fonction pour extraire le profil à une position donnée
function extract_profile(solver, mesh, lx, x_position)
    # Identifier l'indice x le plus proche de x_position
    x_idx = round(Int, (x_position - mesh.x0[1]) / mesh.h[1][1])
    

    # Extraire les valeurs des différentes phases à cette position
    u1ₒ = solver.x[1:(end ÷ 4)]
    u1ᵧ = solver.x[(end ÷ 4 + 1):(end ÷ 2)]
    u2ₒ = solver.x[(end ÷ 2 + 1):(3 * end ÷ 4)]
    u2ᵧ = solver.x[(3 * end ÷ 4 + 1):end]

    # Reshaper les données
    nx = mesh.nx
    ny = mesh.ny
    u1ₒ = reshape(u1ₒ, (nx + 1, ny + 1))'
    u1ᵧ = reshape(u1ᵧ, (nx + 1, ny + 1))'
    u2ₒ = reshape(u2ₒ, (nx + 1, ny + 1))'
    u2ᵧ = reshape(u2ᵧ, (nx + 1, ny + 1))'

    # Extraire le profil le long de y
    y_coords = mesh.y
    profile = (
        y = y_coords,
        u1ₒ = u1ₒ[ix, :],
        u1ᵧ = u1ᵧ[ix, :],
        u2ₒ = u2ₒ[ix, :],
        u2ᵧ = u2ᵧ[ix, :]
    )

    return profile
end

# Paramètres de la simulation
mesh_resolutions = [40, 80, 160, 320]  # Différentes résolutions de maillage
lx, ly = 8.0, 8.0
x0, y0 = 0.0, 0.0
radius = ly / 4
center = (lx / 2 + 0.01, ly / 2 + 0.01)
Δt = 0.01
Tend = 0.3
x_position = lx / 2.01  # Position pour extraire le profil

# Initial condition (peut être ajustée si nécessaire)
u0ₒ1_init = ones((mesh_resolutions[1] + 1) * (mesh_resolutions[1] + 1))
u0ᵧ1_init = zeros((mesh_resolutions[1] + 1) * (mesh_resolutions[1] + 1))
u0ₒ2_init = zeros((mesh_resolutions[1] + 1) * (mesh_resolutions[1] + 1))
u0ᵧ2_init = zeros((mesh_resolutions[1] + 1) * (mesh_resolutions[1] + 1))
u0_init = vcat(u0ₒ1_init, u0ᵧ1_init, u0ₒ2_init, u0ᵧ2_init)

# Stocker les profils pour chaque maillage
profiles = Dict{Int, Any}()

# Boucle sur les différentes résolutions de maillage
for nx in mesh_resolutions
    ny = nx  # Maillage carré

    println("Résolution du maillage : nx = ny = $nx")

    # Résoudre le problème pour le maillage actuel
    mesh, solver = solve_diffusion_unsteady(nx, ny, lx, ly, x0, y0, radius, center, Δt, Tend, u0_init)

    # Extraire le profil
    profile = extract_profile(solver, mesh, lx, x_position)

    # Stocker le profil
    profiles[nx] = profile

end

# Tracer les profils pour différentes résolutions
function plot_convergence_profiles(profiles, lx, x_position)
    fig = Figure()
    ax1 = Axis(fig[1, 1], title = "Convergence des profils - u1ₒ", xlabel = "y", ylabel = "u1ₒ")
    ax2 = Axis(fig[1, 2], title = "Convergence des profils - u1ᵧ", xlabel = "y", ylabel = "u1ᵧ")
    ax3 = Axis(fig[2, 1], title = "Convergence des profils - u2ₒ", xlabel = "y", ylabel = "u2ₒ")
    ax4 = Axis(fig[2, 2], title = "Convergence des profils - u2ᵧ", xlabel = "y", ylabel = "u2ᵧ")

    colors = [:blue, :green, :orange, :red]

    for (i, nx) in enumerate(sort(keys(profiles)))
        profile = profiles[nx]
        y = profile.y
        plot!(ax1, y, profile.u1ₒ, label = "nx=$nx", color = colors[i])
        plot!(ax2, y, profile.u1ᵧ, label = "nx=$nx", color = colors[i])
        plot!(ax3, y, profile.u2ₒ, label = "nx=$nx", color = colors[i])
        plot!(ax4, y, profile.u2ᵧ, label = "nx=$nx", color = colors[i])
    end

    axislegend(ax1, position = :rt)
    axislegend(ax2, position = :rt)
    axislegend(ax3, position = :rt)
    axislegend(ax4, position = :rt)

    display(fig)
end

# Appeler la fonction de tracé
plot_convergence_profiles(profiles, lx, x_position)
