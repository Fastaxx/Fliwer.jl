using Test
using Fliwer

@testset "Capacity.jl" begin
    # Write your tests here.
    body = AutoBody((x,t)->√sum(abs2, x .- 64) - 16)
    fluid = Fluid{2,Float64}(body)

    println(fluid)

    println(fluid.body)
    println(fluid.A)
    println(fluid.B)
    println(fluid.V)
    println(fluid.W)
end