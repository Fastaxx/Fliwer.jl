using Test
using Fliwer

@testset "AutoBody.jl" begin

    # Define the body 1D
    body = Body((x, _=0) -> x^2 - 1.0,(x,_=0)->(x),((-1.,1.),),false)
    @test body.sdf(1.0,0) == 0.

    

    
end