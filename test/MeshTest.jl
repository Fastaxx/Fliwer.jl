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
    @test mesh.faces == (([0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5],),)

    # Test 2D
    nx, ny = 10, 5
    hx, hy = ones(nx), ones(ny)
    x0, y0 = 0., 0.
    mesh = CartesianMesh((hx,hy), (x0,y0))

    @test nC(mesh) == nx*ny
    @test mesh.centers == ([0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10.], [0., 1., 2., 3., 4., 5.],)
    @test mesh.nodes == ([0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5], [0.5, 1.5, 2.5, 3.5, 4.5],)
    @test mesh.faces[1] == ([0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5], [0.0, 1.0, 2.0, 3.0, 4.0, 5.0])
    @test mesh.faces[2] == ([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0], [0.5, 1.5, 2.5, 3.5, 4.5])

    # Test 3D
    nx, ny, nz = 10, 5, 3
    hx, hy, hz = ones(nx), ones(ny), ones(nz)
    x0, y0, z0 = 0., 0., 0.
    mesh = CartesianMesh((hx,hy,hz), (x0,y0,z0))

    @test nC(mesh) == nx*ny*nz
    @test mesh.centers == ([0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10.], [0., 1., 2., 3., 4., 5.], [0., 1., 2., 3.],)
    @test mesh.nodes == ([0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5], [0.5, 1.5, 2.5, 3.5, 4.5], [0.5, 1.5, 2.5],)
    @test mesh.faces[1] == ([0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5], [0.0, 1.0, 2.0, 3.0, 4.0, 5.0], [0.0, 1.0, 2.0, 3.0])
    @test mesh.faces[2] == ([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0], [0.5, 1.5, 2.5, 3.5, 4.5], [0.0, 1.0, 2.0, 3.0])
    @test mesh.faces[3] == ([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0], [0.0, 1.0, 2.0, 3.0, 4.0, 5.0], [0.5, 1.5, 2.5])

    # Test 2D Uniform CartesianMesh
    nx, ny = 10, 5
    lx, ly = 10., 5.
    x0, y0 = 0., 0.
    mesh = CartesianMesh((nx, ny), (lx, ly), (x0, y0))

    @test nC(mesh) == nx*ny
    @test mesh.centers == ([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0], [0.0, 1.0, 2.0, 3.0, 4.0])
    @test mesh.nodes == ([0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5], [0.5, 1.5, 2.5, 3.5, 4.5, 5.5])
    @test mesh.faces[1] == ([0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5], [0.0, 1.0, 2.0, 3.0, 4.0])
    @test mesh.faces[2] == ([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0], [0.5, 1.5, 2.5, 3.5, 4.5]) 

    circle = Body((x, y, _=0) -> (x-5)^2 + (y-2)^2 - 2.0, (x, y, _=0) -> (x, y, _=0), ((-1.0, 1.0), (-1.0, 1.0)), false)

    identify!(mesh, circle)
    
    @test mesh.tag.cut_cells == Tuple{CartesianIndex, Int64}[(CartesianIndex(4, 1), 4), (CartesianIndex(5, 1), 5), (CartesianIndex(6, 1), 6), (CartesianIndex(4, 2), 14), (CartesianIndex(6, 2), 16), (CartesianIndex(4, 3), 24), (CartesianIndex(5, 3), 25), (CartesianIndex(6, 3), 26)]
end