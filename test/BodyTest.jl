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

    

    
end