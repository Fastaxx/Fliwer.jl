using Test
using Fliwer
using SparseArrays

@testset "CapacityTest.jl" begin
    # Cas Test 1D
    nx = 10
    hx = ones(nx)
    x0 = 0.
    body = Body((x, _=0) -> x^2 - 1.0, (x, _=0) -> (x, _=0), ((-1.0, 1.0),), false)
    mesh = CartesianMesh((hx,), (x0,))

    fluid = Fluid(body, mesh)

    # Vérifier que les longueurs des tuples sont correctes
    @test length(fluid.A) == 1
    @test length(fluid.B) == 1
    @test length(fluid.W) == 1

    # Vérifier que la matrice V a la bonne taille
    @test size(fluid.V) == (nx+1, nx+1)

    # Vérifier les types des champs
    @test typeof(fluid.A) == NTuple{1, SparseMatrixCSC{Float64, Int}}
    @test typeof(fluid.B) == NTuple{1, SparseMatrixCSC{Float64, Int}}
    @test typeof(fluid.W) == NTuple{1, SparseMatrixCSC{Float64, Int}}
    @test typeof(fluid.V) == SparseMatrixCSC{Float64, Int}

    # Cas Test 2D
    nx, ny = 10, 5
    hx, hy = ones(nx), ones(ny)
    x0, y0 = 0., 0.
    body = Body((x, y, _=0) -> x^2 + y^2 - 1.0, (x, y, _=0) -> (x, y, _=0), ((-1.0, 1.0), (-1.0, 1.0)), false)
    mesh = CartesianMesh((hx,hy), (x0,y0))

    fluid = Fluid(body, mesh)
    
    # Vérifier que les longueurs des tuples sont correctes
    @test length(fluid.A) == 2
    @test length(fluid.B) == 2
    @test length(fluid.W) == 2

    # Vérifier que la matrice V a la bonne taille
    @test size(fluid.V) == ((nx+1)*(ny+1), (nx+1)*(ny+1))

    # Vérifier les types des champs
    @test typeof(fluid.A) == NTuple{2, SparseMatrixCSC{Float64, Int}}
    @test typeof(fluid.B) == NTuple{2, SparseMatrixCSC{Float64, Int}}
    @test typeof(fluid.W) == NTuple{2, SparseMatrixCSC{Float64, Int}}
    @test typeof(fluid.V) == SparseMatrixCSC{Float64, Int}


    # Cas Test 3D
    nx, ny, nz = 10, 10, 10
    hx, hy, hz = ones(nx), ones(ny), ones(nz)
    x0, y0, z0 = 0., 0., 0.
    body = Body((x, y, z, _=0) -> (x-1.0)^2 + (y-1.0)^2 + (z-1.0)^2 - 1.0 , (x, y, z, _=0) -> (x, y, z, _=0), ((-1.0, 1.0), (-1.0, 1.0), (-1.0, 1.0)), false)
    mesh = CartesianMesh((hx,hy,hz), (x0,y0,z0))

    fluid = Fluid(body, mesh)

    # Vérifier que les longueurs des tuples sont correctes
    @test length(fluid.A) == 3
    @test length(fluid.B) == 3
    @test length(fluid.W) == 3

    # Vérifier que la matrice V a la bonne taille
    @test size(fluid.V) == ((nx+1)*(ny+1)*(nz+1), (nx+1)*(ny+1)*(nz+1))

    # Vérifier les types des champs
    @test typeof(fluid.A) == NTuple{3, SparseMatrixCSC{Float64, Int}}
    @test typeof(fluid.B) == NTuple{3, SparseMatrixCSC{Float64, Int}}
    @test typeof(fluid.W) == NTuple{3, SparseMatrixCSC{Float64, Int}}
    @test typeof(fluid.V) == SparseMatrixCSC{Float64, Int}

end