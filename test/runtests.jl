using Fliwer
using Test

@testset "Fliwer.jl" begin
    # Write your tests here.
    
    include("MeshTest.jl")

    include("BodyTest.jl")
    include("CapacityTest.jl")
end
