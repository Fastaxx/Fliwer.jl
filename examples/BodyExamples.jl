using Fliwer

# Define the body : Circle
n, m = 3*2^6, 2^7
center, radius = m/2, m/8
body = AutoBody((x,t)->√sum(abs2, x .- center) - radius)

# plot with Makie and CairoMakie the signed distance field
x = range(0,stop=n,length=n)
y = range(0,stop=m,length=m)
a = zeros(n,m)
measure_sdf!(a,body)
heatmap(x,y,a,show_axis=false)

