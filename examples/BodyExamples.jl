using Fliwer

# Define the body : Circle
n, m = 3*2^6, 2^7
center, radius = m/2, m/8
body = AutoBody((x,t)->√sum(abs2, x .- center) - radius)

# Measure the body at x=[1,0] and t=0
d,n,V = measure(body, [1,0], 0)

