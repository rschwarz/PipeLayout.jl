import PipeLayout.JunctionLocation: make_nlp
using SCIP # need solver for nonconvex problems :-\

@testset "solve junction location for three terminals (NLP)" begin
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

    @testset "little flow, smallest diameter, Steiner node in center" begin
        inst = Instance(nodes, demand, bounds, diams, ploss_coeff_nice)
        model, x, y, L, l, π = make_nlp(inst, topo, solver)
        status = solve(model, suppress_warnings=true)

        @test status in [:Optimal, :UserLimit]

        xsol, ysol = getvalue(x), getvalue(y)
        for i=1:3 # fixed terminals
            @test_approx_eq_eps xsol[i] nodes[i].x 0.001
            @test_approx_eq_eps ysol[i] nodes[i].y 0.001
        end
        @test_approx_eq_eps xsol[4] 20 0.01
        @test_approx_eq_eps ysol[4] sqrt(3)/6*40 0.01

        Lsol = getvalue(L)
        for i=1:3
            @test_approx_eq_eps Lsol[i] sqrt(3)/3*40 0.01
        end

        lsol = getvalue(l)
        # smallest diameter
        @test_approx_eq_eps sum(lsol[:,1]) 3 0.01
        # no other
        @test_approx_eq_eps sum(lsol[:, 2:end]) 0 0.01
    end

    @testset "more flow, mixed diameter, Steiner node towards source" begin
        inst = Instance(nodes, 20*demand, bounds, diams, ploss_coeff_nice)
        model, x, y, L, l, π = make_nlp(inst, topo, solver)
        status = solve(model, suppress_warnings=true)

        @test status in [:Optimal, :UserLimit]

        xsol, ysol = getvalue(x), getvalue(y)
        for i=1:3 # fixed terminals
            @test_approx_eq_eps xsol[i] nodes[i].x 0.001
            @test_approx_eq_eps ysol[i] nodes[i].y 0.001
        end
        @test_approx_eq_eps xsol[4] 20 0.1
        @test ysol[4] >= sqrt(3)/6*40 # move near source

        Lsol = getvalue(L)
        @test Lsol[1] <= sqrt(3)/3*40
        @test Lsol[2] >= sqrt(3)/3*40
        @test Lsol[3] >= sqrt(3)/3*40

        lsol = getvalue(l)
        D = [d.value for d in diams]
        equiv = (lsol * D.^(-5)).^(-1/5)

        # symmetry
        @test_approx_eq_eps equiv[2] equiv[3] 0.1
        # more flow -> larger diam?
        @test equiv[1] >= equiv[2]
    end
end
