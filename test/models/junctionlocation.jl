using PipeLayout.JunctionLocation
import PipeLayout.JunctionLocation: make_soc
import PipeLayout: ploss_coeff_nice

using JuMP
using SCS

if Pkg.installed("SCIP") != nothing
   include("junctionlocation_nlp.jl")
end

@testset "solve junction location for three terminals (SOC)" begin
    # equilateral triangle
    nodes = [Node(0,0), Node(40,0), Node(20, sqrt(3)/2*40)]
    demand = [10, 10, -20]
    bounds = fill(Bounds(40, 80), 3)
    diams = [Diameter(t...) for t in [(0.6, 0.5), (0.8, 0.8),(1.0, 1.2)]]

    # single Steiner node with dummy location
    topo = Topology(vcat(nodes, [Node(20, 20)]),
                    [Arc(3,4), Arc(4,1), Arc(4,2)])

    solver = SOC(SCSSolver(eps=1e-8, verbose=0))

    @testset "little flow, smallest diameter, Steiner node in center" begin
        inst = Instance(nodes, demand, bounds, diams, ploss_coeff_nice)
        model, x, y, t, π = make_soc(inst, topo, solver)
        status = solve(model, suppress_warnings=true)

        @test status in [:Optimal, :UserLimit]

        xsol, ysol = getvalue(x), getvalue(y)
        for i=1:3 # fixed terminals
            @test_approx_eq_eps xsol[i] nodes[i].x 0.001
            @test_approx_eq_eps ysol[i] nodes[i].y 0.001
        end
        @test_approx_eq_eps xsol[4] 20 0.01
        @test_approx_eq_eps ysol[4] sqrt(3)/6*40 0.01

        tsol = getvalue(t)
        @test_approx_eq_eps sum(tsol[1]) sum(tsol[2]) 0.01
        @test_approx_eq_eps sum(tsol[1]) sum(tsol[3]) 0.01
    end

    @testset "more flow, mixed diameter, Steiner node towards source" begin
        inst = Instance(nodes, 20*demand, bounds, diams, ploss_coeff_nice)
        model, x, y, t, π = make_soc(inst, topo, solver)
        status = solve(model, suppress_warnings=true)

        @test status == :Optimal

        xsol, ysol = getvalue(x), getvalue(y)
        for i=1:3 # fixed terminals
            @test_approx_eq_eps xsol[i] nodes[i].x 0.001
            @test_approx_eq_eps ysol[i] nodes[i].y 0.001
        end
        @test_approx_eq_eps xsol[4] 20 0.1
        @test ysol[4] >= sqrt(3)/6*40 # move near source


        tsol = getvalue(t)
        @test_approx_eq_eps tsol[2] tsol[3] 0.01
    end
end
