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

    # Identify cells
    identify!(mesh, body)

    capacity = Capacity(body, mesh)

    # Vérifier que les longueurs des tuples sont correctes
    @test length(capacity.A) == 1
    @test length(capacity.B) == 1
    @test length(capacity.W) == 1

    # Vérifier que la matrice V a la bonne taille
    @test size(capacity.V) == (nx, nx)

    # Vérifier les types des champs
    @test typeof(capacity.A) == NTuple{1, SparseMatrixCSC{Float64, Int}}
    @test typeof(capacity.B) == NTuple{1, SparseMatrixCSC{Float64, Int}}
    @test typeof(capacity.W) == NTuple{1, SparseMatrixCSC{Float64, Int}}
    @test typeof(capacity.V) == SparseMatrixCSC{Float64, Int}

    # Cas Test 2D
    nx, ny = 10, 5
    hx, hy = ones(nx), ones(ny)
    x0, y0 = 0., 0.
    body = Body((x, y, _=0) -> x^2 + y^2 - 1.0, (x, y, _=0) -> (x, y, _=0), ((-1.0, 1.0), (-1.0, 1.0)), false)
    mesh = CartesianMesh((hx,hy), (x0,y0))

    capacity = Capacity(body, mesh)
    
    # Vérifier que les longueurs des tuples sont correctes
    @test length(capacity.A) == 2
    @test length(capacity.B) == 2
    @test length(capacity.W) == 2

    # Vérifier que la matrice V a la bonne taille
    @test size(capacity.V) == ((nx)*(ny), (nx)*(ny))

    # Vérifier les types des champs
    @test typeof(capacity.A) == NTuple{2, SparseMatrixCSC{Float64, Int}}
    @test typeof(capacity.B) == NTuple{2, SparseMatrixCSC{Float64, Int}}
    @test typeof(capacity.W) == NTuple{2, SparseMatrixCSC{Float64, Int}}
    @test typeof(capacity.V) == SparseMatrixCSC{Float64, Int}


    # Cas Test 3D - Fix CartesianGeometry before running this test
#=     nx, ny, nz = 10, 10, 10
    hx, hy, hz = ones(nx), ones(ny), ones(nz)
    x0, y0, z0 = 0., 0., 0.
    body = Body((x, y, z, _=0) -> (x-1.0)^2 + (y-1.0)^2 + (z-1.0)^2 - 1.0 , (x, y, z, _=0) -> (x, y, z, _=0), ((-1.0, 1.0), (-1.0, 1.0), (-1.0, 1.0)), false)
    mesh = CartesianMesh((hx,hy,hz), (x0,y0,z0))

    capacity = Capacity(body, mesh)

    # Vérifier que les longueurs des tuples sont correctes
    @test length(capacity.A) == 3
    @test length(capacity.B) == 3
    @test length(capacity.W) == 3

    # Vérifier que la matrice V a la bonne taille
    @test size(capacity.V) == ((nx+1)*(ny+1)*(nz+1), (nx+1)*(ny+1)*(nz+1))

    # Vérifier les types des champs
    @test typeof(capacity.A) == NTuple{3, SparseMatrixCSC{Float64, Int}}
    @test typeof(capacity.B) == NTuple{3, SparseMatrixCSC{Float64, Int}}
    @test typeof(capacity.W) == NTuple{3, SparseMatrixCSC{Float64, Int}}
    @test typeof(capacity.V) == SparseMatrixCSC{Float64, Int} =#

end