using Test
using Fliwer
using SparseArrays

@testset "OperatorsTest.jl" begin
    # Test des opérateurs élémentaires
    nx = 5
    periodicity = false

    Dm = ẟ_m(nx, periodicity)
    Dp = δ_p(nx, periodicity)
    Sm = Σ_m(nx, periodicity)
    Sp = Σ_p(nx, periodicity)
    I_op = I(nx)

    @test size(Dm) == (nx, nx)
    @test size(Dp) == (nx, nx)
    @test size(Sm) == (nx, nx)
    @test size(Sp) == (nx, nx)
    @test size(I_op) == (nx, nx)

    # Test de la structure DiffusionOps pour 1D
    A = [spdiagm(0 => ones(nx))]
    B = [spdiagm(0 => ones(nx))]
    V = spdiagm(0 => ones(nx))
    W = (spdiagm(0 => ones(nx)),)

    diffusion_ops_1d = DiffusionOps(A, B, V, W, (nx,))
    @test typeof(diffusion_ops_1d) == DiffusionOps{1}

    # Test de la structure DiffusionOps pour 2D
    nx, ny = 3, 4
    n = nx * ny

    Dx_m = kron(I(ny), ẟ_m(nx))
    Dy_m = kron(ẟ_m(ny), I(nx))

    A = (sprand(n, n, 0.5), sprand(n, n, 0.5))
    B = (sprand(n, n, 0.5), sprand(n, n, 0.5))
    V = sprand(n, n, 0.5)
    W = (sprand(n, n, 0.5), sprand(n, n, 0.5))

    diffusion_ops_2d = DiffusionOps(A, B, V, W, (nx, ny))

    @test typeof(diffusion_ops_2d) == DiffusionOps{2}

    # Test de la structure DiffusionOps pour 3D
    nx, ny, nz = 2, 3, 4
    n = nx * ny * nz

    Dx_m = kron(I(nz), kron(I(ny), ẟ_m(nx)))
    Dy_m = kron(I(nz), kron(ẟ_m(ny), I(nx)))
    Dz_m = kron(ẟ_m(nz), kron(I(ny), I(nx)))

    A = (sprand(n, n, 0.5), sprand(n, n, 0.5), sprand(n, n, 0.5))
    B = (sprand(n, n, 0.5), sprand(n, n, 0.5), sprand(n, n, 0.5))
    V = sprand(n, n, 0.5)
    W = (sprand(n, n, 0.5), sprand(n, n, 0.5), sprand(n, n, 0.5))

    diffusion_ops_3d = DiffusionOps(A, B, V, W, (nx, ny, nz))

    @test typeof(diffusion_ops_3d) == DiffusionOps{3}

    # Test de la fonction ConvectionOps 1D
    A = [spdiagm(0 => ones(nx))]
    B = [spdiagm(0 => ones(nx))]
    V = spdiagm(0 => ones(nx))
    W = (spdiagm(0 => ones(nx)))

    uₒ, uᵧ = rand(nx), rand(nx)

#    convection_ops_1d = ConvectionOps(A, B, V, W, (nx,), uₒ, uᵧ)

    #@test typeof(convection_ops_1d) == ConvectionOps{1}

    # Test de la fonction ConvectionOps 2D
    nx, ny = 3, 4
    n = nx * ny
    A = (sprand(n, n, 0.5), sprand(n, n, 0.5))
    B = (sprand(n, n, 0.5), sprand(n, n, 0.5))
    V = sprand(n, n, 0.5)
    W = (sprand(n, n, 0.5), sprand(n, n, 0.5))
    
    uₒx, uₒy = rand(n), rand(n)
    uᵧx, uᵧy = rand(n), rand(n)

    uₒ, uᵧ = vcat(uₒx, uₒy), vcat(uᵧx, uᵧy)

#    convection_ops_2d = ConvectionOps(A, B, V, W, (nx, ny), uₒ, uᵧ)

  #  @test typeof(convection_ops_2d) == ConvectionOps{2}

    # Test de la fonction ConvectionOps 3D
    nx, ny, nz = 2, 3, 4
    n = nx * ny * nz
    A = (sprand(n, n, 0.5), sprand(n, n, 0.5), sprand(n, n, 0.5))
    B = (sprand(n, n, 0.5), sprand(n, n, 0.5), sprand(n, n, 0.5))
    V = sprand(n, n, 0.5)
    W = (sprand(n, n, 0.5), sprand(n, n, 0.5), sprand(n, n, 0.5))

    uₒx, uₒy, uₒz = rand(n), rand(n), rand(n)
    uᵧx, uᵧy, uᵧz = rand(n), rand(n), rand(n)

    uₒ, uᵧ = vcat(uₒx, uₒy, uₒz), vcat(uᵧx, uᵧy, uᵧz)

   # convection_ops_3d = ConvectionOps(A, B, V, W, (nx, ny, nz), uₒ, uᵧ)

  #  @test typeof(convection_ops_3d) == ConvectionOps{3}
   
end