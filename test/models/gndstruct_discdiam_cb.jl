@testset "Solve semi decomposition with nogoods on y (callback)" begin
    #       7    9      even arc numbers for
    #   () - d2 - ()    reversed arcs
    #   /1   /3   /5
    #  s1 - () - d1
    #    11   13

    # slightly irregular grid
    nodes = [Node(0,0), Node(45,0), Node(25, 22)]
    demand = [-50, 20, 30]
    bounds = fill(Bounds(60,80), 3)
    diams = [Diameter(t...) for t in [(0.8, 1.0),(1.0, 1.2)]]

    topo = Topology([Node(t...) for t in [(0,22), (0,0), (25,22), (25,0), (45,22), (45,0)]],
                    [Arc(t...) for t in [(1,2), (2,1), (3,4), (4,3), (5,6), (6,5),
                                         (1,3), (3,1), (3,5), (5,3),
                                         (2,4), (4,2), (4,6), (6,4)]])

    mastersolver = SCIP.Optimizer(display_verblevel=0)
    subsolver = SCIP.Optimizer(display_verblevel=0)
    cbtoposolver = GndStr.CallbackTopo(mastersolver, subsolver)

    @testset "low flow: very easy instance" begin
        inst = Instance(nodes, 1 * demand, bounds, diams)
        result = optimize(inst, topo, cbtoposolver)
        @test result.status == :Optimal

        zsol = result.solution.zsol
        qsol = result.solution.qsol

        # shortest tree is obvious, smallest diameter enough:
        @test zsol[11, 1] ≈ 1.0
        @test zsol[13, 1] ≈ 1.0
        @test zsol[4, 1] ≈ 1.0
        @test sum(zsol) ≈ 3.0 # all others 0

        # uniq flow solution
        @test qsol[11] ≈ 50.0
        @test qsol[13] ≈ 20.0
        @test qsol[4] ≈ 30.0
        @test sum(qsol) ≈ 100.0 # all others 0
    end

    @testset "medium flow: difficult instance" begin
        inst = Instance(nodes, 10 * demand, bounds, diams)
        result = optimize(inst, topo, cbtoposolver)
        @test result.status == :Optimal

        zsol = result.solution.zsol
        qsol = result.solution.qsol

        # shortest tree is obvious, need one large diameter
        @test sum(zsol[2, :]) == 1
        @test sum(zsol[7, :]) == 1
        @test sum(zsol[11, :]) == 1
        @test sum(zsol[13, :]) == 1
        @test sum(zsol[:, 1]) == 1
        @test sum(zsol[:, 2]) == 3
        @test sum(zsol) ≈ 4.0 # all others 0
    end

    @testset "high flow on triangle: infeasible instance" begin
        inst = Instance([Node(0,0), Node(50,0)],
                        [-1000, 1000],
                        [Bounds(60,80), Bounds(60,80)],
                        [Diameter(t...) for t in [(0.8, 1.0),(1.0, 1.2)]])
        topo = Topology([Node(0,0), Node(50,0), Node(30, 40)],
                        [Arc(1,3), Arc(1,2), Arc(2,3),
                         Arc(3,1), Arc(2,1), Arc(3,2)])

        result = optimize(inst, topo, cbtoposolver)
        @test result.status == :Infeasible
    end
end


@testset "Solve GBD (callback)" begin
    #       7    9      even arc numbers for
    #   () - d2 - ()    reversed arcs
    #   /1   /3   /5
    #  s1 - () - d1
    #    11   13

    # slightly irregular grid
    nodes = [Node(0,0), Node(45,0), Node(25, 22)]
    demand = [-50, 20, 30]
    bounds = fill(Bounds(60,80), 3)
    diams = [Diameter(t...) for t in [(0.8, 1.0),(1.0, 1.2)]]

    topo = Topology([Node(t...) for t in [(0,22), (0,0), (25,22), (25,0), (45,22), (45,0)]],
                    [Arc(t...) for t in [(1,2), (2,1), (3,4), (4,3), (5,6), (6,5),
                                         (1,3), (3,1), (3,5), (5,3),
                                         (2,4), (4,2), (4,6), (6,4)]])

    mastersolver = SCIP.Optimizer(display_verblevel=0)
    subsolver = GLPK.Optimizer()
    cbgbdsolver = GndStr.CallbackGBD(mastersolver, subsolver)

    @testset "low flow: very easy instance" begin
        inst = Instance(nodes, 1 * demand, bounds, diams)
        result = optimize(inst, topo, cbgbdsolver)
        @test result.status == :Optimal

        zsol = result.solution.zsol
        qsol = result.solution.qsol

        # shortest tree is obvious, smallest diameter enough:
        @test zsol[11, 1] ≈ 1.0
        @test zsol[13, 1] ≈ 1.0
        @test zsol[4, 1] ≈ 1.0
        @test sum(zsol) ≈ 3.0 # all others 0

        # uniq flow solution
        @test qsol[11] ≈ 50.0
        @test qsol[13] ≈ 20.0
        @test qsol[4] ≈ 30.0
        @test sum(qsol) ≈ 100.0 # all others 0
    end

    @testset "medium flow: difficult instance" begin
        inst = Instance(nodes, 10 * demand, bounds, diams)
        result = optimize(inst, topo, cbgbdsolver)
        @test result.status == :Optimal

        zsol = result.solution.zsol
        qsol = result.solution.qsol

        # shortest tree is obvious, need one large diameter
        @test sum(zsol[2, :]) == 1
        @test sum(zsol[7, :]) == 1
        @test sum(zsol[11, :]) == 1
        @test sum(zsol[13, :]) == 1
        @test sum(zsol[:, 1]) == 1
        @test sum(zsol[:, 2]) == 3
        @test sum(zsol) ≈ 4.0 # all others 0
    end

    @testset "high flow on triangle: infeasible instance" begin
        inst = Instance([Node(0,0), Node(50,0)],
                        [-1000, 1000],
                        [Bounds(60,80), Bounds(60,80)],
                        [Diameter(t...) for t in [(0.8, 1.0),(1.0, 1.2)]])
        topo = Topology([Node(0,0), Node(50,0), Node(30, 40)],
                        [Arc(1,3), Arc(1,2), Arc(2,3),
                         Arc(3,1), Arc(2,1), Arc(3,2)])

        result = optimize(inst, topo, cbgbdsolver)
        @test result.status == :Infeasible
    end
end
