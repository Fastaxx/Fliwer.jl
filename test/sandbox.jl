using Fliwer
using Test

# Test Moving Body
@testset "BodyTest.jl" begin
# Test 1D
nx = 10
hx = ones(nx)
x0 = 0.
lx = 10.
mesh = CartesianMesh((hx,), (x0,))
domain = ((x0, lx),)

x_pos = 5.5
map = (x, t)->(x + 0.1 * t)
body = Body((x,_=0)->-(x - x_pos), map, domain, true)

@test sdf(body, 5.5, 0) == 0.0

println("SDF: ", body.sdf(5.5, 0))

d,n,v = measure(body, [4.3,], 1)

println("Distance: ", d)
println("Normal: ", n)
println("Velocity: ", v)

# Test 2D
nx, ny = 10, 5
hx, hy = ones(nx), ones(ny)
x0, y0 = 0., 0.
lx, ly = 5., 5.
domain = ((x0, lx), (y0, ly))
mesh = CartesianMesh((hx, hy), (x0, y0))

radius, center = ly/4, (lx/2, ly/2)
mapping = (t, x, y) -> (x + 0.1 * t, y)
circle = Body((x,y,_=0)->-(sqrt((x-center[1])^2 + (y-center[2])^2) - radius), mapping, domain, true)

@test sdf(circle, 2.5, 2.5) == 1.25

d,n,v = measure(circle, [2.5, 2.5], 1)

println("Distance: ", d)
println("Normal: ", n)
println("Velocity: ", v)


end











"""

# Test Mesh 
@testset "MeshTest.jl" begin
# Test 1D
nx = 10
hx = ones(nx)
x0 = 0.
mesh = CartesianMesh((hx,), (x0,))

@test mesh.faces == (([0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5],),)

# Test 1D : Uniform CartesianMesh
nx = 10
lx = 10.
x0 = 0.
mesh = CartesianMesh((nx,), (lx,), (x0,))

@test mesh.faces == (([0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5],),)

# Test 2D
nx, ny = 10, 5
hx, hy = ones(nx), ones(ny)
x0, y0 = 0., 0.
mesh = CartesianMesh((hx, hy), (x0, y0))

@test mesh.faces[1] == ([0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5], [0.0, 1.0, 2.0, 3.0, 4.0, 5.0])
@test mesh.faces[2] == ([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0], [0.5, 1.5, 2.5, 3.5, 4.5])

@show mesh.centers
@show mesh.nodes
@show mesh.faces

# Test 2D : Uniform CartesianMesh

nx, ny = 10, 5
lx, ly = 10., 5.
x0, y0 = 0., 0.
mesh = CartesianMesh((nx, ny), (lx, ly), (x0, y0))

@show mesh.centers
@show mesh.nodes
@show mesh.faces

@test mesh.faces[1] == ([0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5], [0.0, 1.0, 2.0, 3.0, 4.0, 5.0])
@test mesh.faces[2] == ([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0], [0.5, 1.5, 2.5, 3.5, 4.5])

# Test 3D
nx, ny, nz = 10, 5, 3
hx, hy, hz = ones(nx), ones(ny), ones(nz)
x0, y0, z0 = 0., 0., 0.
mesh = CartesianMesh((hx, hy, hz), (x0, y0, z0))

@test mesh.faces[1] == ([0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5], [0.0, 1.0, 2.0, 3.0, 4.0, 5.0], [0.0, 1.0, 2.0, 3.0])
@test mesh.faces[2] == ([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0], [0.5, 1.5, 2.5, 3.5, 4.5], [0.0, 1.0, 2.0, 3.0])
@test mesh.faces[3] == ([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0], [0.0, 1.0, 2.0, 3.0, 4.0, 5.0], [0.5, 1.5, 2.5])

end
"""