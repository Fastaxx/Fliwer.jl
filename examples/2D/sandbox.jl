using Fliwer
using IterativeSolvers
using SparseArrays
using LinearAlgebra
using CairoMakie

### 2D Test Case : Monophasic Steady Diffusion Equation inside a Disk
# Define the mesh
nx, ny = 80, 80
lx, ly = 4., 4.
x0, y0 = 0., 0.
domain = ((x0, lx), (y0, ly))
mesh = CartesianMesh((nx, ny), (lx, ly), (x0, y0))

# Define the body
radius, center = ly/4, (lx/2, ly/2) .+ (0.01, 0.01)
circle = Body((x,y,_=0)->(sqrt((x-center[1])^2 + (y-center[2])^2) - radius), (x,y,_)->(x,y), domain, false)

# Identify cells
identify!(mesh, circle)

# Define the capacity
capacity = Capacity(circle, mesh)

# Define the operators
operator = DiffusionOps(capacity.A, capacity.B, capacity.V, capacity.W, (nx+1, ny+1))

Ax = capacity.A[1]
Ay = capacity.A[2]
Dx_m = kron(spdiagm(0 => ones(nx+1)), ẟ_m(ny+1))
Dy_m = kron(ẟ_m(ny+1), spdiagm(0 => ones(nx+1)))

Iγp = sqrt.((Ax*Dx_m).^2 + (Ay*Dy_m).^2)
Iγp = diag(Iγp)
# Remove the max value
Iγp[Iγp .> maximum(Iγp)*0.99999999999999] .= 0.0

function build_I_g(operator::AbstractOperators)
    n = prod(operator.size)
    vec_1 = [1 for i in 1:size(operator.H', 2)]
    row_sums = operator.H'*vec_1
    abs_row_sums = abs.(row_sums)
    #abs_row_sums[abs_row_sums .< 10^-8] .= 0.0  # Mettre à zéro les petites valeurs
    Iᵧ = spdiagm(0 => abs_row_sums)
    return Iᵧ
end

Iᵧ = build_I_g(operator)
Γ = diag(capacity.Γ)

# Reshape 
Iᵧ = diag(Iᵧ)
Iᵧ = reshape(Iᵧ, (nx+1, ny+1))
Γ = reshape(Γ, (nx+1, ny+1))
Iγp = reshape(Iγp, (nx+1, ny+1))


# Plot 
using CairoMakie

fig = Figure(size = (800, 800))
ax = Axis(fig[1, 1], aspect = DataAspect(), xlabel = "x", ylabel = "y", title="Iᵧ")
ax1 = Axis(fig[1, 2], aspect = DataAspect(), xlabel = "x", ylabel = "y", title="Γ")
ax2 = Axis(fig[1, 3], aspect = DataAspect(), xlabel = "x", ylabel = "y", title="Iγp")
hm = heatmap!(ax, Iᵧ, colormap = :viridis, label="Iᵧ")
hm2 = heatmap!(ax1, Γ, colormap = :viridis, label="Γ")
hm3 = heatmap!(ax2, Iγp, colormap = :viridis, label="Iγp")
Colorbar(fig[1, 4], hm)
Colorbar(fig[1, 5], hm2)
Colorbar(fig[1, 6], hm3)
fig
