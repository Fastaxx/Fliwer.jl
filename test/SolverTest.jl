using Test
using Fliwer

@testset "SolverTest.jl" begin
    # Write your tests here.

    s1 = DiffusionSolver(true, true, 1.0)
    initialize(s1)

    @test s1.time_type == TimeType(0)
    @test s1.phase_type == PhaseType(0)
    @test s1.equation_type == EquationType(0)

end
