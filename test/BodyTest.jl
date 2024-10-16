using Test
using Fliwer

@testset "AutoBody Tests" begin
    sdf(x, t) = sum(x) + t
    map(x, t) = x .+ t

    @testset "AutoBody constructor" begin
        body = AutoBody(sdf, map)
        @test body.sdf([1, 2], 3) == 6
        @test body.map([1, 2], 3) == [4, 5]
    end

    @testset "AutoBody operations" begin
        body1 = AutoBody(sdf, map)
        body2 = AutoBody((x, t) -> sum(x) - t, map)

        @testset "AutoBody addition" begin
            body3 = body1 + body2
            @test body3.sdf([1, 2], 3) == 0
            @test body3.map([1, 2], 3) == [4, 5]
        end

        @testset "AutoBody intersection" begin
            body4 = body1 ∩ body2
            @test body4.sdf([1, 2], 3) == 6
            @test body4.map([1, 2], 3) == [4, 5]
        end
    end
end

@testset "Bodies Tests" begin
    sdf(x, t) = sum(x) + t
    map(x, t) = x .+ t
    body1 = AutoBody(sdf, map)
    body2 = AutoBody((x, t) -> sum(x) - t, map)

    @testset "Bodies constructor" begin
        bodies = Bodies([body1, body2], [Base.:+])
        @test bodies.bodies[1].sdf([1, 2], 3) == 6
        @test bodies.bodies[2].sdf([1, 2], 3) == 0
    end

    @testset "Bodies operations" begin
        bodies1 = Bodies([body1, body2], [Base.:+])
        bodies2 = Bodies([body1, body2], [Base.:-])

        @testset "Bodies addition" begin
            bodies3 = bodies1 + bodies2
            @test length(bodies3.bodies) == 4
            @test length(bodies3.ops) == 3
        end
    end
end