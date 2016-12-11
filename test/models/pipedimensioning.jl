using PipeLayout.PipeDimensioning
using JuMP
using Clp

@testset "check dimensioning on single pipe" begin
    inst, topo = single_pipe()
    solver = PipeDimLP(ClpSolver())

    result = optimize(inst, topo, solver)
    @test result.status == :Optimal

    πsol = result.sol.πsol
    psol = πsol.^0.5
    @test approx_in(psol[1], 60.0, 80.0)
    @test approx_in(psol[2], 60.0, 80.0)
    @test psol[1] >= psol[2]

    lsol = result.sol.lsol
    @test size(lsol) == (1,3) # one pipe with 3 segments
    @test sum(lsol[1,:]) ≈ 1.0
    @test approx_in(lsol[1,1], 0.0, 1.0)
    @test approx_in(lsol[1,2], 0.0, 1.0)
    @test approx_in(lsol[1,3], 0.0, 1.0)
    @test is_SOS2(lsol[1,:])
end
