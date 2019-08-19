using PipeLayout.PipeDim
import PipeLayout: isapproxin, isSOS2
using JuMP
using GLPK

@testset "check dimensioning on single pipe" begin
    inst, topo = single_pipe()
    solver = PipeDim.LP(GLPK.Optimizer)

    result = optimize(inst, topo, solver)
    @test result.status == :Optimal

    πsol = result.sol.πsol
    psol = πsol.^0.5
    @test isapproxin(psol[1], 60.0, 80.0)
    @test isapproxin(psol[2], 60.0, 80.0)
    @test psol[1] > psol[2]

    lsol = result.sol.lsol
    @test size(lsol) == (1,3) # one pipe with 3 segments
    @test sum(lsol[1,:]) ≈ 1.0
    @test isapproxin(lsol[1,1], 0.0, 1.0)
    @test isapproxin(lsol[1,2], 0.0, 1.0)
    @test isapproxin(lsol[1,3], 0.0, 1.0)
    @test isSOS2(lsol[1,:])
end
