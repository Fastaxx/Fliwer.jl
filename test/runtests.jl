using Fliwer
using Test

@testset "Fliwer.jl" begin
    # Write your tests here.
    
    include("MeshTest.jl")
    include("BodyTest.jl")
    include("CapacityTest.jl")
    include("OperatorsTest.jl")
    include("UtilsTest.jl")
    include("BoundaryTest.jl")
    #include("SolverTest.jl")
end
