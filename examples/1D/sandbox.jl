using Fliwer
using IterativeSolvers
using SparseArrays, StaticArrays
using LinearAlgebra, SpecialFunctions
using CairoMakie
using DelimitedFiles, Roots

# Analytical solution for the Stefan problem
T₀ = 1.0    # Hot temperature
Tm = 0.0    # Cold temperature
L = 1.0     # Latent heat
c = 1.0     # Thermal capacity
k = 1.0     # Thermal conductivity
ρ = 1.0     # Density

Stefan_number = c*T₀/L
println("Stefan number: ", Stefan_number)

f = (λ) -> λ*exp(λ^2)*erf(λ) - Stefan_number/sqrt(π)

λ = find_zero(f, 0.0)
println("λ: ", λ)
function T(x, t)
    return T₀ - T₀/erf(λ) * erf(x/(2*sqrt(k*t)))
end

function ∇T(x,t)
    return -T₀/(sqrt(t)) * exp(-x^2/(4*t)) / (sqrt(π) * erf(λ))
end

function s(t)
    return 2*λ*sqrt(k*t)
end

function v(t)
    return λ*sqrt(k/t)
end

# Plot of $T(x,t)$ at times $t=[0.25,0.5,0.75,1]$
x = range(0, stop=1.0, length=200)
t = [0.25, 0.5, 0.75, 1.0]
y = [T.(x, tᵢ) for tᵢ in t]

y[1][x .>= s(t[1])] .= 0.0
y[2][x .>= s(t[2])] .= 0.0
y[3][x .>= s(t[3])] .= 0.0
y[4][x .>= s(t[4])] .= 0.0

# Plot the temperature
fig = Figure()
ax = Axis(fig[1, 1]; xlabel="x", ylabel="T(x,t)", title="Stefan problem - Temperature - St=$Stefan_number")
lines!(ax, x, y[1], color=:blue, linewidth=2, label="t=0.25")
lines!(ax, x, y[2], color=:red, linewidth=2, label="t=0.5")
lines!(ax, x, y[3], color=:green, linewidth=2, label="t=0.75")
lines!(ax, x, y[4], color=:black, linewidth=2, label="t=1.0")
axislegend(ax)
display(fig)
readline()
# Plot the gradient of the temperature
y = [∇T.(x, tᵢ) for tᵢ in t]

fig = Figure()
ax = Axis(fig[1, 1]; xlabel="x", ylabel="∇T(x,t)", title="Stefan problem - Gradient - St=$Stefan_number")
lines!(ax, x, y[1], color=:blue, linewidth=2, label="t=0.25")
lines!(ax, x, y[2], color=:red, linewidth=2, label="t=0.5")
lines!(ax, x, y[3], color=:green, linewidth=2, label="t=0.75")
lines!(ax, x, y[4], color=:black, linewidth=2, label="t=1.0")
axislegend(ax, position=:rb)
display(fig)

readline()
# Plot the interface position
t=range(0, stop=1.0, length=200)
fig = Figure()
ax = Axis(fig[1, 1]; xlabel="t", ylabel="s(t)", title="Stefan problem - Interface position - St=$Stefan_number")
lines!(ax, t, s.(t), color=:blue, linewidth=2, label="Interface position")
axislegend(ax, position=:rb)
display(fig)
readline()

# Plot the velocity of the interface
fig = Figure()
ax = Axis(fig[1, 1]; xlabel="t", ylabel="v(t)", title="Stefan problem - Interface velocity - St=$Stefan_number")
lines!(ax, t, v.(t), color=:blue, linewidth=2, label="Interface velocity")
axislegend(ax)
display(fig)