import PipeLayout.PipeDimensioning: make_model
using JuMP
using Clp

@testset "check dimensioning on single pipe" begin
    inst, topo = single_pipe()
    model, π, l = make_model(inst, topo, ClpSolver())
    status = solve(model)
    @test status == :Optimal

    psol = getvalue(π).^0.5
    @test approx_in(psol[1], 60.0, 80.0)
    @test approx_in(psol[2], 60.0, 80.0)
    @test psol[1] >= psol[2]

    lsol = getvalue(l)
    @test size(lsol) == (1,3) # one pipe with 3 segments
    @test sum(lsol[1,:]) ≈ 1.0
    @test approx_in(lsol[1,1], 0.0, 1.0)
    @test approx_in(lsol[1,2], 0.0, 1.0)
    @test approx_in(lsol[1,3], 0.0, 1.0)
    @test is_SOS2(lsol[1,:])
end
