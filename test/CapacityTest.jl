using Test
using Fliwer
using SparseArrays

@testset "Test Fluid Construction" begin
    # Cas Test 1D
    nx = 10
    hx = ones(nx)
    x0 = 0.
    body = AutoBody((x,t)->√sum(abs2, x .- 0.5) - 0.25)
    mesh = CartesianMesh((hx,), (x0,))

    fluid = Fluid(body, mesh)

    println("Fluid : ", fluid)

    # Vérifier que les longueurs des tuples sont correctes
    @test length(fluid.A) == 1
    @test length(fluid.B) == 1
    @test length(fluid.W) == 1

    # Vérifier que la matrice V a la bonne taille
    @test size(fluid.V) == (nx, nx)

    # Vérifier les types des champs
    @test typeof(fluid.A) == NTuple{1, SparseMatrixCSC{Float64, Int}}
    @test typeof(fluid.B) == NTuple{1, SparseMatrixCSC{Float64, Int}}
    @test typeof(fluid.W) == NTuple{1, SparseMatrixCSC{Float64, Int}}
    @test typeof(fluid.V) == SparseMatrixCSC{Float64, Int}
end