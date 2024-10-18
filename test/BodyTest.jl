using Test
using Fliwer

@testset "AutoBody.jl" begin
    norm2(x) = √sum(abs2,x)

    # test AutoDiff in 2D and 3D
    body1 = AutoBody((x,t)->norm2(x)-2-t)
    @test all(measure(body1,[√2.,√2.],0.).≈(0,[√.5,√.5],[0.,0.]))
    @test all(measure(body1,[2.,0.,0.],1.).≈(-1.,[1.,0.,0.],[0.,0.,0.]))

    body2 = AutoBody((x,t)->norm2(x)-2,(x,t)->x.+t^2)
    @test all(measure(body2,[√2.,√2.],0.).≈(0,[√.5,√.5],[0.,0.]))
    @test all(measure(body2,[1.,-1.,-1.],1.).≈(0.,[1.,0.,0.],[-2.,-2.,-2.]))

    
end