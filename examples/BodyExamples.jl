using Fliwer

# Define the body : Circle
n, m = 3*2^6, 2^7
center, radius = m/2, m/8
body = AutoBody((x,t)->√sum(abs2, x .- center) - radius)

# Measure the body at x=[1,0] and t=0
d,no,V = measure(body, [1,0], 0)

# Print the results
println("d = $d")
println("no = $no")
println("V = $V")

# Plot the body
a = zeros(n,m)
for i ∈ 1:n, j ∈ 1:m
    a[i,j] = measure(body, [j,i], 0)[1]
end
println("a = $a")

using Plots
Plots.default(show=true)
Plots.heatmap(a', aspect_ratio=1, c=:viridis)
readline()