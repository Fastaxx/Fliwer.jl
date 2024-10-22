using Test
using Fliwer

@testset "MeshTest.jl" begin
    # Write your tests here.
    
    # Test 1D
    nx = 10
    hx = ones(nx)
    x0 = 0.
    mesh = CartesianMesh((hx,), (x0,))

    @test nC(mesh) == nx
    @test mesh.centers == ([0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10.],)
    @test mesh.nodes == ([0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5],)

    # Test 2D
    nx, ny = 10, 5
    hx, hy = ones(nx), ones(ny)
    x0, y0 = 0., 0.
    mesh = CartesianMesh((hx,hy), (x0,y0))

    @test nC(mesh) == nx*ny
    @test mesh.centers == ([0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10.], [0., 1., 2., 3., 4., 5.],)
    @test mesh.nodes == ([0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5], [0.5, 1.5, 2.5, 3.5, 4.5],)

    # Test 3D
    nx, ny, nz = 10, 5, 3
    hx, hy, hz = ones(nx), ones(ny), ones(nz)
    x0, y0, z0 = 0., 0., 0.
    mesh = CartesianMesh((hx,hy,hz), (x0,y0,z0))

    @test nC(mesh) == nx*ny*nz
    @test mesh.centers == ([0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10.], [0., 1., 2., 3., 4., 5.], [0., 1., 2., 3.],)
    @test mesh.nodes == ([0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5], [0.5, 1.5, 2.5, 3.5, 4.5], [0.5, 1.5, 2.5],)

end