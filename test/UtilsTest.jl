using Test
using Fliwer

@testset "UtilsTest.jl" begin
    # Test find_border
    # Test 1D
    mesh_1d = CartesianMesh((ones(10),), (0.,))
    border_cells_1d = find_border(mesh_1d)

    @test length(border_cells_1d) == 2
    @test border_cells_1d[1] == (CartesianIndex(1), 1)
    @test border_cells_1d[2] == (CartesianIndex(11), 11)

    # Test 2D
    mesh_2d = CartesianMesh((ones(10), ones(5)), (0., 0.))
    border_cells_2d = find_border(mesh_2d)

    @test length(border_cells_2d) == 30
    @test border_cells_2d[1] == (CartesianIndex(1, 1), 1)
    @test border_cells_2d[30] == (CartesianIndex(11, 6), 66)
    @test border_cells_2d[15] == (CartesianIndex(11, 3), 33)

    # Test 3D
    mesh_3d = CartesianMesh((ones(10), ones(5), ones(3)), (0., 0., 0.))
    border_cells_3d = find_border(mesh_3d)

    @test length(border_cells_3d) == 192
    @test border_cells_3d[1] == (CartesianIndex(1, 1, 1), 1)
    @test border_cells_3d[60] == (CartesianIndex(5, 6, 1), 60)
    @test border_cells_3d[30] == (CartesianIndex(8, 3, 1), 30)


    # Test find_cut
    # Test 1D
    body = Body((x, _=0) -> x - 5.0,(x,_=0)->(x),((-1.,1.),),false)
    sdf_1D = eval_sdf(mesh_1d, body)
    cut_cells_1d = find_cut(mesh_1d, body)

    @test length(cut_cells_1d) == 1
    @test cut_cells_1d[1] == (CartesianIndex(5,), 5)

    # Test 2D
    body = Body((x, y, _=0) -> (x-5)^2 + (y-2)^2 - 1.0, (x, y, _=0) -> (x, y, _=0), ((-1.0, 1.0), (-1.0, 1.0)), false)
    sdf_2D = eval_sdf(mesh_2d, body)
    cut_cells_2d = find_cut(mesh_2d, body)

    @test length(cut_cells_2d) == 8
    @test cut_cells_2d[1] == (CartesianIndex(4, 1), 4)

    tag = MeshTag([], [], [])
    identify!(tag, mesh_2d, body)

end