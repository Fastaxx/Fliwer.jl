using Test
using Fliwer

@testset "AutoBody.jl" begin

    # Define the body 1D
    body = Body((x, _=0) -> x^2 - 1.0,(x,_=0)->(x),((-1.,1.),),false)
    @test body.sdf(1.0,0) == 0.

    # Define the body 2D 
    body1 = Body((x,y,_=0)->sqrt(x^2 + y^2) - 1, (x,y,_)->(x,y), ((-1.0, 1.0), (-1.0, 1.0)), false)
    body2 = Body((x,y,_=0)->sqrt(x^2 + y^2) - 0.5, (x,y,_)->(x,y), ((-1.0, 1.0), (-1.0, 1.0)), false)

    body = body1 + body2
    @test sdf(body1, 0.0, 0.0) == -1.0
    @test sdf(body2, 0.0, 0.0) == -0.5
    @test sdf(body, 0.0, 0.0) == -1.0

    body = body1 ⊖ body2
    @test sdf(body, 0.0, 0.0) == 0.5

    body = c(body1)
    @test sdf(body, 0.0, 0.0) == 1.0

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