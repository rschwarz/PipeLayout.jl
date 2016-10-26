using PipeLayout.JunctionLocation
import PipeLayout: ploss_coeff_nice
import PipeLayout.JunctionLocation: make_nlp

using JuMP
using SCIP # need solver for nonconvex problems :-\

facts("solve junction location for three terminals") do
    # equilateral triangle
    nodes = [Node(0,0), Node(40,0), Node(20, sqrt(3)*20)]
    demand = [10, 10, -20]
    bounds = fill(Bounds(40, 80), 3)
    diams = [Diameter(t...) for t in [(0.6, 0.5), (0.8, 0.8),(1.0, 1.2)]]

    # single Steiner node with dummy location
    topo = Topology(vcat(nodes, [Node(20, 20)]),
                    [Arc(3,4), Arc(4,1), Arc(4,2)])

    solver = NLP(SCIPSolver("display/verblevel", 0))

    context("little flow, smallest diameter, node in center") do
        inst = Instance(nodes, demand, bounds, diams, ploss_coeff_nice)
        model, x, y, L, l, π = make_nlp(inst, topo, solver)
        status = solve(model)

        @fact status --> :Optimal

        lsol = getvalue(l)
        @fact sum(lsol) --> roughly(3.0)
        @fact sum(lsol[:,1]) --> roughly(3.0) # smallest diameter
    end
end
