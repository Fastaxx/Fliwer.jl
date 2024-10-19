using Fliwer
using Test

@testset "Capacity.jl" begin
    # Write your tests here.
    body = NoBody()
    fluid = Fluid{2,Float64}(body)

    println(fluid)

    println(fluid.body)
    println(fluid.A)
    println(fluid.B)
    println(fluid.V)
    println(fluid.W)
end