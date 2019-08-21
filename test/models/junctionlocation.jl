using PipeLayout.JuncLoc

using JuMP
using SCS

@testset "solve junction location for three terminals (SOC)" begin
    # equilateral triangle
    nodes = [Node(0,0), Node(40,0), Node(20, sqrt(3)/2*40)]
    demand = [10, 10, -20]
    bounds = fill(Bounds(40, 80), 3)
    diams = [Diameter(t...) for t in [(0.6, 0.5), (0.8, 0.8),(1.0, 1.2)]]

    # single Steiner node with dummy location
    topo = Topology(vcat(nodes, [Node(20, 20)]),
                    [Arc(3,4), Arc(4,1), Arc(4,2)])

    solver = JuncLoc.SOC(SCSSolver(eps=1e-6, verbose=0))

    @testset "little flow, smallest diameter, Steiner node in center" begin
        inst = Instance(nodes, demand, bounds, diams, ploss_coeff_nice)
        model, x, y, t, π = JuncLoc.make_soc(inst, topo, solver)
        JuMP.optimize!(model)
	    status = JuMP.termination_status(model)

        @test status in [:Optimal, :UserLimit]

        xsol, ysol = getvalue(x), getvalue(y)
        for i=1:3 # fixed terminals
            @test xsol[i] ≈ nodes[i].x atol=0.001
            @test ysol[i] ≈ nodes[i].y atol=0.001
        end
        @test xsol[4] ≈ 20 atol=0.01
        @test ysol[4] ≈ sqrt(3)/6*40 atol=0.01

        tsol = getvalue(t)
        @test sum(tsol[1]) ≈ sum(tsol[2]) atol=0.01
        @test sum(tsol[1]) ≈ sum(tsol[3]) atol=0.01
    end

    @testset "more flow, mixed diameter, Steiner node towards source" begin
        inst = Instance(nodes, 20*demand, bounds, diams, ploss_coeff_nice)
        model, x, y, t, π = JuncLoc.make_soc(inst, topo, solver)
        JuMP.optimize!(model)
	    status = JuMP.termination_status(model)

        @test status == :Optimal

        xsol, ysol = getvalue(x), getvalue(y)
        for i=1:3 # fixed terminals
            @test xsol[i] ≈ nodes[i].x atol=0.001
            @test ysol[i] ≈ nodes[i].y atol=0.001
        end
        @test xsol[4] ≈ 20 atol=0.1
        @test ysol[4] >= sqrt(3)/6*40 # move near source


        tsol = getvalue(t)
        @test tsol[2] ≈ tsol[3] atol=0.01
    end

    @testset "using the optimize function to solve" begin
        inst = Instance(nodes, 20*demand, bounds, diams, ploss_coeff_nice)
        result = optimize(inst, topo, solver)
        @test result.status == :Optimal

        sol = result.sol
        toposol = Topology(sol.nodes, topo.arcs)
        L = pipelengths(toposol)
        c = [d.cost for d in diams]
        obj = L' * sol.lsol * c
        @test result.value ≈ obj[1] atol=0.001
    end
end
