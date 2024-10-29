function plot_solution(u, nx, ny)
    reshaped_u = reshape(u[1:length(u) ÷ 2], (nx + 1, ny + 1))'
    fig = Figure()
    ax = Axis(fig[1, 1], title = "Solution Plot", xlabel = "x", ylabel = "y")
    hm = heatmap!(ax, reshaped_u, colormap = :viridis)
    Colorbar(fig[1, 2], hm, label = "Intensity")
    display(fig)
end
