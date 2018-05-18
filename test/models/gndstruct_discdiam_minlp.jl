@testset "run MINLP model" begin
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

    @testset "low flow: very easy instance" begin
        inst = Instance(nodes, 1*demand, bounds, diams)

        result = optimize(inst, topo, GndStr.MINLP(SCIPSolver("display/verblevel", 0)))
        @test result.status == :Optimal

        zsol = result.solution.zsol
        @test zsol[4,1] == true
        @test zsol[11,1] == true
        @test zsol[13,1] == true
        @test sum(zsol) == 3

        @test result.dualbound ≈ 67.0
    end

    @testset "medium flow: difficult instance" begin
        inst = Instance(nodes, 5*demand, bounds, diams)

        result = optimize(inst, topo, GndStr.MINLP(SCIPSolver("display/verblevel", 0)))
        @test result.status == :Optimal

        zsol = result.solution.zsol
        @test zsol[4,1] == true
        @test zsol[11,2] == true
        @test zsol[13,1] == true
        @test sum(zsol) == 3

        @test result.dualbound ≈ 72.0
    end

    @testset "high flow on triangle: infeasible" begin
        inst3 = Instance([Node(0,0), Node(50,0)],
                         20*[-50, 50],
                         [Bounds(60,80), Bounds(60,80)],
                         [Diameter(t...) for t in [(0.8, 1.0),(1.0, 1.2)]])
        topo3 = Topology([Node(0,0), Node(50,0), Node(30, 40)],
                         [Arc(1,3), Arc(1,2), Arc(2,3)])

        result = optimize(inst3, topo3, GndStr.MINLP(SCIPSolver("display/verblevel", 0)))
        @test result.status == :Infeasible
        @test result.solution == nothing
        @test result.dualbound == Inf
    end

    @testset "on the example with subproblem/relaxation gap" begin
        #     __s1__  __s2
        #   t3      t4
        nodes = [Node(100,0), Node(300,0), Node(0,0), Node(200,0)]
        arcs = [Arc(1,3), Arc(1,4), Arc(2,4)]
        topo = Topology(nodes, arcs)

        demand = [-600, -400, 400, 600]
        bounds = fill(Bounds(40, 80), 4)
        diams = [Diameter(1.0, 1.0), Diameter(2.0, 2.0)]
        inst = Instance(nodes, demand, bounds, diams, ploss_coeff_nice)

        result = optimize(inst, topo, GndStr.MINLP(SCIPSolver("display/verblevel", 0)))
        @test result.status == :Optimal
        zsol = result.solution.zsol
        @test sum(zsol[:,2]) == 1
        @test sum(zsol[:,1]) == 2 # solution not quite uniqe
        @test result.dualbound ≈ 400.0
    end

    @testset "difficult instance with disconnected candidates" begin
        inst = Instance([Node(0,0), Node(100,0), Node(200,100)],
                        [800, -900, 100],
                        fill(Bounds(40,80), 3),
                        [Diameter(t...) for t in [(1.0, 1.0), (2.0, 3.2)]],
                        ploss_coeff_nice)
        topo = squaregrid(2, 3, 100.0, antiparallel=true)

        result = optimize(inst, topo, GndStr.MINLP(SCIPSolver("display/verblevel", 0)))
        @test result.status == :Optimal

        zsol = result.solution.zsol
        @test sum(zsol) == 3

        @test result.dualbound ≈ 520.0
    end
end
