using Test
using Fliwer

@testset "BoundaryTest.jl" begin
    nx, ny = 40, 40
    lx, ly = 4.0, 4.0
    x0, y0 = 0.0, 0.0
    mesh = CartesianMesh((nx, ny), (lx, ly), (x0, y0))

    radius, center = ly/4, (lx/2, ly/2)
    circle = Body((x,y,_=0)->sqrt((x-center[1])^2 + (y-center[2])^2) - radius, (x,y,_)->(x,y), ((x0, lx), (y0, ly)), false)

    identify!(mesh, circle)

    borders = Dict{Symbol, AbstractBoundary}()

    borders[:left] = Dirichlet(0.0)
    borders[:right] = Dirichlet(0.0)
    borders[:top] = Dirichlet(0.0)
    borders[:bottom] = Dirichlet(0.0)

    interfaces = Dict{Symbol, AbstractInterfaceBC}()

    scalar = ScalarJump(1.0, 1.0, 0.0)
    flux = FluxJump(1.0, 1.0, 0.0)

    bc = BorderConditions(borders)
    ic = InterfaceConditions(scalar, flux)
    
    @test ic == InterfaceConditions(ScalarJump(1.0, 1.0, 0.0), FluxJump(1.0, 1.0, 0.0))

end