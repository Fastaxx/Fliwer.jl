using Test
using Fliwer
using SparseArrays

@testset "SolverTest.jl" begin
    # Write your tests here.
    nx, ny = 40, 40
    lx, ly = 4.0, 4.0
    x0, y0 = 0.0, 0.0
    domain=((x0, lx), (y0, ly))
    mesh = CartesianMesh((nx, ny), (lx, ly), (x0, y0))  # Exemple 2D

    # Define the body
    radius, center = ly/4, (lx/2, ly/2) .+ (0.01, 0.01)
    circle = Body((x,y,_=0)->sqrt((x-center[1])^2 + (y-center[2])^2) - radius, (x,y,_)->(x,y), domain, false)

    # Identify cells
    identify!(mesh, circle)

    bc = Dirichlet(1.0)
    bc_b = BorderConditions(Dict{Symbol, AbstractBoundary}(:left => bc, :right => bc, :top => bc, :bottom => bc))

    A=spzeros((nx+1)*(ny+1), (nx+1)*(ny+1))
    b=zeros((nx+1)*(ny+1))

    

end
