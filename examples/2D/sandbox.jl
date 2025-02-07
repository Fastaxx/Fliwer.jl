using Fliwer
using IterativeSolvers
using SparseArrays, StaticArrays
using LinearAlgebra, SpecialFunctions
using CairoMakie

using CairoMakie
using DelimitedFiles

# Analytical Solution : 2d Heat Equation inside a disk 2ph Henry
# Analytical solution
using QuadGK
using SpecialFunctions

Dg,Dl = 1.0, 1.0
R0 = 1.0
cg0, cl0 = 1.0, 0.0
He = 1.0

D = sqrt(Dg/Dl)

function Phi(u)
    term1 = Dg*sqrt(Dl)*besselj1(u*R0)*bessely0(D*u*R0)
    term2 = He*Dl*sqrt(Dg)*besselj0(u*R0)*bessely1(D*u*R0)
    return term1 - term2
end

function Psi(u)
    term1 = Dg*sqrt(Dl)*besselj1(u*R0)*besselj0(D*u*R0)
    term2 = He*Dl*sqrt(Dg)*besselj0(u*R0)*besselj1(D*u*R0)
    return term1 - term2
end

function cg_integrand(u, r, t)
    Φu = Phi(u)
    Ψu = Psi(u)
    denom = u^2*(Φu^2 + Ψu^2)
    num   = exp(-Dg*u^2*t)*besselj0(u*r)*besselj1(u*R0)
    return iszero(denom) ? 0.0 : num/denom
end

function cl_integrand(u, r, t)
    Φu = Phi(u)
    Ψu = Psi(u)
    denom = u*(Φu^2 + Ψu^2)
    term1 = besselj0(D*u*r)*Φu
    term2 = bessely0(D*u*r)*Ψu
    num   = exp(-Dg*u^2*t)*besselj1(u*R0)*(term1 - term2)
    return iszero(denom) ? 0.0 : num/denom
end

function compute_cg(r_values, t_values)
    prefactor = (4*cg0*Dg*Dl*Dl*He)/(π^2*R0)
    cg_results = Matrix{Float64}(undef, length(t_values), length(r_values))
    for (i, t) in pairs(t_values)
        Umax = 5.0/sqrt(Dg*t)
        for (j, r) in pairs(r_values)
            val, _ = quadgk(u->cg_integrand(u, r, t), 0, Umax; atol=1e-6, rtol=1e-6)
            cg_results[i, j] = prefactor*val
        end
    end
    return cg_results
end

function compute_cl(r_values, t_values)
    prefactor = (2*cg0*Dg*sqrt(Dl)*He)/π
    cl_results = Matrix{Float64}(undef, length(t_values), length(r_values))
    for (i, t) in pairs(t_values)
        Umax = 5.0/sqrt(Dg*t)
        for (j, r) in pairs(r_values)
            val, _ = quadgk(u->cl_integrand(u, r, t), 0, Umax; atol=1e-6, rtol=1e-6)
            cl_results[i, j] = prefactor*val
        end
    end
    return cl_results
end

function cg_once(r, t)
    prefactor = (4*cg0*Dg*Dl*Dl*He)/(π^2*R0)
    Umax = 5.0 / sqrt(Dg * t)
    val, _ = quadgk(u -> cg_integrand(u, r, t), 0, Umax; atol=1e-6, rtol=1e-6)
    return prefactor * val
end

function cl_once(r, t)
    prefactor = (2*cg0*Dg*sqrt(Dl)*He)/π
    Umax = 5.0 / sqrt(Dg * t)
    val, _ = quadgk(u -> cl_integrand(u, r, t), 0, Umax; atol=1e-6, rtol=1e-6)
    return prefactor * val
end

xvals = range(-2, 2, length=200)
yvals = range(-2, 2, length=200)
t     = 1.0999999999999999  # pick a single time for visualization

cg_2d = [cg_once(sqrt(x^2 + y^2), t) for y in yvals, x in xvals]
cl_2d = [cl_once(sqrt(x^2 + y^2), t) for y in yvals, x in xvals]

# After computing cg_2d and cl_2d, do:
for (j, y) in enumerate(yvals)
    for (i, x) in enumerate(xvals)
        r = sqrt(x^2 + y^2)
        if r > R0
            # Outside
            cg_2d[j, i] = NaN
        else
            # Inside
            cl_2d[j, i] = NaN
        end
    end
end

# Plot the analytical solution
using CairoMakie

fig = Figure()
ax1 = Axis(fig[1,1], aspect=DataAspect(), xlabel="x", ylabel="y", title="Phase 1: cg(x,y)")
hm1 = heatmap!(ax1, xvals, yvals, cg_2d, colormap=:viridis)
Colorbar(fig[1,2], hm1, label="cg")

ax2 = Axis(fig[2,1], aspect=DataAspect(), xlabel="x", ylabel="y", title="Phase 2: cl(x,y)")
hm2 = heatmap!(ax2, xvals, yvals, cl_2d, colormap=:plasma)
Colorbar(fig[2,2], hm2, label="cl")

display(fig)

readline()

# Poisson Equation inside a square
x0, y0 = 0.0, 0.0
lx, ly = 4.0, 4.0
nx, ny = 40, 40
x=range(x0, stop=lx, length=nx+1)
y=range(y0, stop=ly, length=ny+1)

radius, center = ly/4.0, (lx/2, ly/2)

# Analytical solution
function u_analytical(x,y)
    # Sum from m,n=1 to 100 m,n odd of 16 L²/((π^4 mn) * (m² + n²)) * sin(mπx/L) * sin(nπy/L)
    x = x #- center[1]
    y = y #- center[2]
    sum = 0.0
    for m in 1:2:500
        for n in 1:2:500
            sum += 16 * lx^2 / (π^4 * m * n * (m^2 + n^2)) * sin(m*π*x/lx) * sin(n*π*y/ly)
        end
    end
    return sum
end

u_ana = [u_analytical(x_, y_) for x_ in x, y_ in y]

fig = Figure()
ax = Axis(fig[1, 1], aspect = DataAspect(), xlabel = "x", ylabel = "y", title="Analytical Solution")
hm = heatmap!(ax, u_ana, colormap = :viridis)
Colorbar(fig[1, 2], hm)
display(fig)

println(maximum(u_ana))

readline()


# Divergence operator
# Define the mesh
nx, ny = 40, 40
lx, ly = 4., 4.
x0, y0 = 0., 0.
domain = ((x0, lx), (y0, ly))
mesh = CartesianMesh((nx, ny), (lx, ly), (x0, y0))

# Define the body
radius, center = ly/4.0, (lx/2, ly/2) #.+ (0.01, 0.01)
circle = Body((x,y,_=0)->(sqrt((x-center[1])^2 + (y-center[2])^2) - radius), (x,y,_)->(x,y), domain, false)

# Identify cells
identify!(mesh, circle)

# Define the capacity
capacity = Capacity(circle, mesh)

# Define the operators
operator = DiffusionOps(capacity.A, capacity.B, capacity.V, capacity.W, (nx+1, ny+1))

# Build Divergence Operators
x=range(x0, stop=lx, length=nx+1)
y=range(y0, stop=ly, length=ny+1)

x_faces = (x[1:end-1] .+ x[2:end]) ./ 2
y_faces = (y[1:end-1] .+ y[2:end]) ./ 2

Coord_facesx = [(x_faces[i], y[j]) for i in 1:length(x_faces), j in 1:length(y)]
Coord_facesy = [(x[i], y_faces[j]) for i in 1:length(x), j in 1:length(y_faces)]

qxω = [capacity.C_ω[i][1] for i in 1:length(capacity.C_ω)]
qyω = [capacity.C_ω[i][2] for i in 1:length(capacity.C_ω)]
qxγ = [capacity.C_γ[i][1] for i in 1:length(capacity.C_γ)]
qyγ = [capacity.C_γ[i][2] for i in 1:length(capacity.C_γ)]
qω = vcat(qxω, qyω)
qγ = vcat(qxγ, qyγ)

Divergence = ∇_(operator, qω, qγ)
Divergence = [Divergence[i] * capacity.V[i,i] for i in 1:length(Divergence)]
Divergence = reshape(Divergence, (nx+1, ny+1))
Divergence[capacity.cell_types .== 0] .= NaN

println(Divergence)

using CairoMakie

fig = Figure()
ax = Axis(fig[1, 1], aspect = DataAspect(), xlabel = "x", ylabel = "y", title="Divergence")
hm = heatmap!(ax, Divergence, colormap = :viridis)
Colorbar(fig[1, 2], hm)
display(fig)

readline()

# Volume Redefinition
# Define the mesh
nx = 20
lx = 4.0
x0 = 0.0
domain = ((x0, lx),)
mesh = CartesianMesh((nx,), (lx,), (x0,))

# Define the body
center, radius = lx/2, 0.5
body = Body((x,_=0)->sqrt((x - center)^2) - radius, (x,_)->(x), domain, false)

# Identify cells
identify!(mesh, body)

# Define the capacity
capacity = Capacity(body, mesh)

# Define the operators
operator = DiffusionOps(capacity.A, capacity.B, capacity.V, capacity.W, (nx+1,))

# Redefine W volume
xₒ = capacity.C_ω[1]
pₒ = [capacity.C_ω[i][1] for i in 1:length(capacity.C_ω)]
pᵧ = [capacity.C_γ[i][1] for i in 1:length(capacity.C_ω)]
p = vcat(pₒ, pᵧ)
grad = ∇(operator, p)
W_new = [grad[i] * capacity.W[1][i,i] for i in 1:length(grad)]
W_new = spdiagm(0 => W_new)
@show W_new
@show capacity.W[1]
# In 1D the new W volume is the same as the old one
readline()

# Redefine V volume : Quadratic Profile in the gradient then divergence with the gradient
xₒ = capacity.C_ω[1]
pₒ = [(capacity.C_ω[i][1]^2)/2 for i in 1:length(capacity.C_ω)]
pᵧ = [(capacity.C_γ[i][1]^2)/2 for i in 1:length(capacity.C_ω)]

p = vcat(pₒ, pᵧ)
grad = ∇(operator, p)

qω = vcat(grad)
qγ = vcat(grad)

div = ∇_(operator, qω, qγ)
V_new = spdiagm(0 => div)
@show V_new
@show capacity.V

"""
fig = Figure()
ax = Axis(fig[1, 1], xlabel = "x", ylabel = "Grad")
scatter!(ax, grad, color=:blue, label="Grad")
scatter!(ax, grad_true, color=:red, label="Grad True")
axislegend(ax, position=:rb)
display(fig)

readline()

error = abs.(grad - grad_true)./grad_true
error[capacity.cell_types .== 0] .= NaN
fig=Figure()
ax = Axis(fig[1, 1], xlabel = "x", ylabel = "Error")
scatter!(ax, error, color=:blue, label="Error")
axislegend(ax, position=:rb)
display(fig)
"""

readline()


# Moving 2D Analytical Solution

# λ is the root of the equation λ² exp(λ²) Ei(-λ²) +(Tₘ - Tₒ)/L = 0 : Function
Tₘ = 0.0
Tₒ = 1.0
L = 1.0
D = 1.0
R₀ = 1.0
λ = λ_0(Tₘ, Tₒ, L)
f = (λ) -> λ^2 * exp(λ^2) * expinti(-λ^2) + (Tₘ - Tₒ) / L

function λ_0(Tₘ, Tₒ, L)
    f = (λ) -> λ^2 * exp(λ^2) * expinti(-λ^2) + (Tₘ - Tₒ) / L
    return find_zero(f, 0.0)
end

# Plot the function
λs = range(-1000.0, stop=1000.0, length=1000)
fs = [f(λ) for λ in λs]

fig = Figure()
ax = Axis(fig[1, 1], xlabel = "λ", ylabel = "f(λ)", title="Root of the equation")
lines!(ax, λs, fs, color=:blue)
display(fig)
readline()
    
R(t) = R₀ + 2λ*√(D*t)

T(x, y, t) = Tₒ + ((Tₘ - Tₒ) * expinti(-(x^2 + y^2)/(4*D)*t))/(expinti(-λ^2))

# Plot the solution
x = range(0, stop=2, length=100)
y = range(0, stop=2, length=100)
t = 1.0

T_ = [T(x_, y_, t) for x_ in x, y_ in y]

fig = Figure()
ax = Axis(fig[1, 1], aspect = DataAspect(), xlabel = "x", ylabel = "y", title="Analytical Solution")
hm = heatmap!(ax, T_, colormap = :viridis)
Colorbar(fig[1, 2], hm)
display(fig)


readline()
# Interface Centroid
# Define the mesh
nx, ny = 20, 20
lx, ly = 4., 4.
x0, y0 = 0., 0.
domain = ((x0, lx), (y0, ly))
mesh = CartesianMesh((nx, ny), (lx, ly), (x0, y0))

# Define the body
radius, center = ly/4.0, (lx/2, ly/2) #.+ (0.01, 0.01)
circle = Body((x,y,_=0)->(sqrt((x-center[1])^2 + (y-center[2])^2) - radius), (x,y,_)->(x,y), domain, false)

# Identify cells
identify!(mesh, circle)

# Define the capacity
capacity = Capacity(circle, mesh)

@show size(capacity.C_ω)
@show size(capacity.C_γ)

x_coords,y_coords = mesh.nodes
Φ = circle.sdf
Φ = (x) -> (x[1] - center[1])^2 + (x[2] - center[2])^2 - radius^2

f = circle.sdf
# Convert f(x,y,_=0) to a function of the form Φ(x) with x = (x,y)
f = (x) -> circle.sdf(x[1], x[2], 0.0)

"""
# Interface centroids
using ImplicitIntegration

function computeInterfaceCentroids(mesh::CartesianMesh{1},body)
    x_coords = mesh.nodes[1]
    nx = length(x_coords)-1
    Φ = (x) -> body.sdf(x[1], 0.0, 0.0)


    C_γ = [0.0 for _ in 1:nx]
    Γ   = zeros(nx)

    for i in 1:nx
        a = x_coords[i]
        b = x_coords[i+1]

        area_val = ImplicitIntegration.integrate((p)->1, Φ, a, b, surface=true).val
        if area_val > 0
            x_c = ImplicitIntegration.integrate((p)->p[1], Φ, a, b, surface=true).val / area_val
            C_γ[i] = x_c
            Γ[i]   = area_val
        else
            C_γ[i] = 0.0
            Γ[i]   = 0.0
        end
    end
    return C_γ, Γ
end

function computeInterfaceCentroids(mesh::CartesianMesh{2}, body)
    x_coords, y_coords = mesh.nodes
    nx, ny = length(x_coords)-1, length(y_coords)-1
    Φ = (x) -> body.sdf(x[1], x[2], 0.0)


    C_γ = [SVector{2, Float64}(0.0, 0.0) for _ in 1:nx, _ in 1:ny]
    Γ   = zeros(nx, ny)

    for i in 1:nx
        for j in 1:ny
            a = (x_coords[i],   y_coords[j])
            b = (x_coords[i+1], y_coords[j+1])

            area_val = ImplicitIntegration.integrate((p)->1, Φ, a, b, surface=true).val
            if area_val > 0
                x_c = ImplicitIntegration.integrate((p)->p[1], Φ, a, b, surface=true).val / area_val
                y_c = ImplicitIntegration.integrate((p)->p[2], Φ, a, b, surface=true).val / area_val
                C_γ[i, j] = SVector(x_c, y_c)
                Γ[i, j]   = area_val
            else
                C_γ[i, j] = SVector(0.0, 0.0)
                Γ[i, j]   = 0.0
            end
        end
    end
    return C_γ, Γ
end

function computeInterfaceCentroids(mesh::CartesianMesh{3}, body)
    x_coords, y_coords, z_coords = mesh.nodes
    nx, ny, nz = length(x_coords)-1, length(y_coords)-1, length(z_coords)-1
    Φ = (x) -> body.sdf(x[1], x[2], x[3])

    C_γ = [SVector{3, Float64}(0.0, 0.0, 0.0) for _ in 1:nx, _ in 1:ny, _ in 1:nz]
    Γ   = zeros(nx, ny, nz)

    for i in 1:nx
        for j in 1:ny
            for k in 1:nz
                a = (x_coords[i],   y_coords[j],   z_coords[k])
                b = (x_coords[i+1], y_coords[j+1], z_coords[k+1])

                area_val = ImplicitIntegration.integrate((p)->1, Φ, a, b, surface=true).val
                if area_val > 0
                    x_c = ImplicitIntegration.integrate((p)->p[1], Φ, a, b, surface=true).val / area_val
                    y_c = ImplicitIntegration.integrate((p)->p[2], Φ, a, b, surface=true).val / area_val
                    z_c = ImplicitIntegration.integrate((p)->p[3], Φ, a, b, surface=true).val / area_val
                    C_γ[i, j, k] = SVector(x_c, y_c, z_c)
                    Γ[i, j, k]   = area_val
                else
                    C_γ[i, j, k] = SVector(0.0, 0.0, 0.0)
                    Γ[i, j, k]   = 0.0
                end
            end
        end
    end
    return C_γ, Γ
end
"""
#C_γ, Γ = computeInterfaceCentroids(mesh, circle)

@show size(capacity.C_γ)
@show size(capacity.Γ)


readline()


# Suppose we have 4 files, each representing a different mesh size.
# For example:
files = [
    "//home/libat/github/Fliwer.jl/max_T_log30.txt"
    "/home/libat/github/Fliwer.jl/max_T_log60.txt"
    "/home/libat/github/Fliwer.jl/max_T_log80.txt"
    "/home/libat/github/Fliwer.jl/max_T_log100.txt"
]

mesh_labels = ["Mesh 30", "Mesh 60", "Mesh 80", "Mesh 100"]
colors = [:blue, :red, :green, :orange]

fig = Figure()
ax = Axis(fig[1, 1], xlabel="Time Step", ylabel="Max(T)", title="Mesh Influence on Oscillations")
max_norm = []
L2_err = []
for (i, file) in enumerate(files)
    data = readdlm(file) |> vec  # Each row is a single Float64 value
    push!(max_norm, maximum(data .- 1.0))
    push!(L2_err, norm(data .- 1.0,2))
    steps = 1:length(data)
    lines!(ax, steps, data, color=colors[i], label=mesh_labels[i])
end

axislegend(ax, position=:rb)
display(fig)

fig = Figure()
ax = Axis(fig[1, 1], xlabel="Mesh Size", ylabel="Max(T)", title="Mesh Influence on Oscillations")
scatter!(ax, [30, 60, 80, 100], max_norm, color=:blue, label="Max(T)")
scatter!(ax, [30, 60, 80, 100], L2_err, color=:red, label="L2 error")
axislegend(ax, position=:rb)
display(fig)


readline()


path_cn = "/home/libat/Bureau/CutCell/TestsCases/TranslatingDisk/max_T_log_CN.txt"
path_be = "/home/libat/Bureau/CutCell/TestsCases/TranslatingDisk/max_T_log_BE.txt"

# Plot the data with CairoMakie
using CairoMakie

path_cn = "/home/libat/Bureau/CutCell/TestsCases/TranslatingDisk/max_T_log_CN.txt"
path_be = "/home/libat/Bureau/CutCell/TestsCases/TranslatingDisk/max_T_log_BE.txt"

# Read data from each file
data_cn = open(path_cn, "r") do io
    lines = readlines(io)
    [parse(Float64, line) for line in lines]
end

data_be = open(path_be, "r") do io
    lines = readlines(io)
    [parse(Float64, line) for line in lines]
end

# Plot on the same figure
fig = Figure()
ax = Axis(fig[1, 1], xlabel = "Index", ylabel = "max(|T|)", title="CN vs BE")

lines!(ax, 1:length(data_cn), data_cn, color=:blue, label="CN")
lines!(ax, 1:length(data_be), data_be, color=:red, label="BE")

axislegend(ax, position=:rb)
display(fig)

readline()









### 2D Test Case : Monophasic Steady Diffusion Equation inside a Disk
# Define the mesh
nx, ny = 80, 80
lx, ly = 4., 4.
x0, y0 = 0., 0.
domain = ((x0, lx), (y0, ly))
mesh = CartesianMesh((nx, ny), (lx, ly), (x0, y0))

# Define the body
radius, center = ly/4.0, (lx/2, ly/2) #.+ (0.01, 0.01)
circle = Body((x,y,_=0)->(sqrt((x-center[1])^2 + (y-center[2])^2) - radius), (x,y,_)->(x,y), domain, false)

# Identify cells
identify!(mesh, circle)

# Define the capacity
capacity = Capacity(circle, mesh)

# Define the operators
operator = DiffusionOps(capacity.A, capacity.B, capacity.V, capacity.W, (nx+1, ny+1))

# Gradient operator
xₒ = capacity.C_ω[1]
pₒ = [capacity.C_ω[i][1]^2 + capacity.C_ω[i][2]^2 for i in 1:length(capacity.C_ω)]
pᵧ = [capacity.C_ω[i][1]^2 + capacity.C_ω[i][2]^2 for i in 1:length(capacity.C_ω)]

p = vcat(pₒ, pᵧ)
∇_ = ∇(operator, p)

@show size(∇_)

∇x = ∇_[1:((nx+1)*(ny+1))]
∇y = ∇_[((nx+1)*(ny+1)+1):end]
∇x[capacity.cell_types .== 0] .= NaN
#∇x[capacity.cell_types .== -1] .= NaN
∇y[capacity.cell_types .== 0] .= NaN
#∇y[capacity.cell_types .== -1] .= NaN


∇x_true = [2*capacity.C_ω[i][1]  for i in 1:length(capacity.C_ω)]
∇y_true = [2*capacity.C_ω[i][2] for i in 1:length(capacity.C_ω)]
∇x_true[capacity.cell_types .== 0] .= NaN
∇y_true[capacity.cell_types .== 0] .= NaN

error_∇x = log10.(abs.((∇x - ∇x_true)./∇x_true))
error_∇y = log10.(abs.((∇y - ∇y_true)./∇y_true))

using CairoMakie

fig = Figure()
ax = Axis(fig[1, 1], aspect = DataAspect(), xlabel = "x", ylabel = "y", title="∇x")
ax1 = Axis(fig[1, 2], aspect = DataAspect(), xlabel = "x", ylabel = "y", title="∇y")
hm = heatmap!(ax, reshape(∇x, (nx+1, ny+1)), colormap = :viridis)
hm2 = heatmap!(ax1, reshape(∇y, (nx+1, ny+1)), colormap = :viridis)
Colorbar(fig[1, 3], hm)
display(fig)
readline()

fig = Figure()
ax = Axis(fig[1, 1], aspect = DataAspect(), xlabel = "x", ylabel = "y", title="Error ∇x")
ax1 = Axis(fig[1, 2], aspect = DataAspect(), xlabel = "x", ylabel = "y", title="Error ∇y")
hm = heatmap!(ax, reshape(error_∇x, (nx+1, ny+1)), colormap = :viridis)
hm2 = heatmap!(ax1, reshape(error_∇y, (nx+1, ny+1)), colormap = :viridis)
Colorbar(fig[1, 3], hm)
display(fig)
readline()



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
