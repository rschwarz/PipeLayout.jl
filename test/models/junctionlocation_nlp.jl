import PipeLayout.JunctionLocation: make_nlp
using SCIP # need solver for nonconvex problems :-\

facts("solve junction location for three terminals (NLP)") do
    # equilateral triangle
    nodes = [Node(0,0), Node(40,0), Node(20, sqrt(3)/2*40)]
    demand = [10, 10, -20]
    bounds = fill(Bounds(40, 80), 3)
    diams = [Diameter(t...) for t in [(0.6, 0.5), (0.8, 0.8),(1.0, 1.2)]]

    # single Steiner node with dummy location
    topo = Topology(vcat(nodes, [Node(20, 20)]),
                    [Arc(3,4), Arc(4,1), Arc(4,2)])

    # don't solve nonconvex NLP to optimality
    solver = NLP(SCIPSolver("display/verblevel", 0,
                            "limits/gap", 1e-3))

    context("little flow, smallest diameter, Steiner node in center") do
        inst = Instance(nodes, demand, bounds, diams, ploss_coeff_nice)
        model, x, y, L, l, π = make_nlp(inst, topo, solver)
        status = solve(model, suppress_warnings=true)

        @fact status --> anyof(:Optimal, :UserLimit)

        xsol, ysol = getvalue(x), getvalue(y)
        for i=1:3 # fixed terminals
            @fact xsol[i] --> roughly(nodes[i].x, 0.001)
            @fact ysol[i] --> roughly(nodes[i].y, 0.001)
        end
        @fact xsol[4] --> roughly(20, 0.01)
        @fact ysol[4] --> roughly(sqrt(3)/6*40, 0.01)

        Lsol = getvalue(L)
        for i=1:3
            @fact Lsol[i] --> roughly(sqrt(3)/3*40, 0.01)
        end

        lsol = getvalue(l)
        @fact sum(lsol[:,1]) --> roughly(3, 0.01) # smallest diameter
        @fact sum(lsol[:, 2:end]) --> roughly(0, 0.01) # no other
    end

    context("more flow, mixed diameter, Steiner node towards source") do
        inst = Instance(nodes, 20*demand, bounds, diams, ploss_coeff_nice)
        model, x, y, L, l, π = make_nlp(inst, topo, solver)
        status = solve(model, suppress_warnings=true)

        @fact status --> anyof(:Optimal, :UserLimit)

        xsol, ysol = getvalue(x), getvalue(y)
        for i=1:3 # fixed terminals
            @fact xsol[i] --> roughly(nodes[i].x, 0.001)
            @fact ysol[i] --> roughly(nodes[i].y, 0.001)
        end
        @fact xsol[4] --> roughly(20, 0.1)
        @fact ysol[4] --> greater_than(sqrt(3)/6*40) # move near source

        Lsol = getvalue(L)
        @fact Lsol[1] --> less_than(sqrt(3)/3*40)
        @fact Lsol[2] --> greater_than(sqrt(3)/3*40)
        @fact Lsol[3] --> greater_than(sqrt(3)/3*40)

        lsol = getvalue(l)
        D = [d.value for d in diams]
        equiv = (lsol * D.^(-5)).^(-1/5)

        @fact equiv[2] --> roughly(equiv[3], 0.1) # symmetry
        @fact equiv[1] --> greater_than(equiv[2]) # more flow -> larger diam?
    end
end
