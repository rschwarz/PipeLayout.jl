using PipeLayout.JunctionLocation
import PipeLayout.JunctionLocation: make_soc
import PipeLayout: ploss_coeff_nice

using JuMP
using SCS

if Pkg.installed("SCIP") != nothing
   include("junctionlocation_nlp.jl")
end

facts("solve junction location for three terminals (SOC)") do
    # equilateral triangle
    nodes = [Node(0,0), Node(40,0), Node(20, sqrt(3)/2*40)]
    demand = [10, 10, -20]
    bounds = fill(Bounds(40, 80), 3)
    diams = [Diameter(t...) for t in [(0.6, 0.5), (0.8, 0.8),(1.0, 1.2)]]

    # single Steiner node with dummy location
    topo = Topology(vcat(nodes, [Node(20, 20)]),
                    [Arc(3,4), Arc(4,1), Arc(4,2)])

    solver = SOC(SCSSolver(eps=1e-8, verbose=0))

    context("little flow, smallest diameter, Steiner node in center") do
        inst = Instance(nodes, demand, bounds, diams, ploss_coeff_nice)
        model, x, y, t, π = make_soc(inst, topo, solver)
        status = solve(model, suppress_warnings=true)

        @fact status --> anyof(:Optimal, :UserLimit)

        xsol, ysol = getvalue(x), getvalue(y)
        for i=1:3 # fixed terminals
            @fact xsol[i] --> roughly(nodes[i].x, 0.001)
            @fact ysol[i] --> roughly(nodes[i].y, 0.001)
        end
        @fact xsol[4] --> roughly(20, 0.01)
        @fact ysol[4] --> roughly(sqrt(3)/6*40, 0.01)

        tsol = getvalue(t)
        @fact sum(tsol[1]) --> roughly(sum(tsol[2]), 0.01)
        @fact sum(tsol[1]) --> roughly(sum(tsol[3]), 0.01)
    end

    context("more flow, mixed diameter, Steiner node towards source") do
        inst = Instance(nodes, 20*demand, bounds, diams, ploss_coeff_nice)
        model, x, y, t, π = make_soc(inst, topo, solver)
        status = solve(model, suppress_warnings=true)

        @fact status --> :Optimal

        xsol, ysol = getvalue(x), getvalue(y)
        for i=1:3 # fixed terminals
            @fact xsol[i] --> roughly(nodes[i].x, 0.001)
            @fact ysol[i] --> roughly(nodes[i].y, 0.001)
        end
        @fact xsol[4] --> roughly(20, 0.1)
        @fact ysol[4] --> greater_than(sqrt(3)/6*40) # move near source


        tsol = getvalue(t)
        @fact tsol[2] --> roughly(tsol[3], 0.01)
    end
end
