@testset "solve master problem: ground structure, discrete diameters" begin
    #       7    9      even arc numbers for
    #   () - d2 - ()    reversed arcs
    #   /1   /3   /5
    #  s1 - () - d1
    #    11   13
    inst = Instance([Node(0,0), Node(40,0), Node(20, 20)],
                    [-50, 20, 30],
                    fill(Bounds(60,80), 3),
                    [Diameter(t...) for t in [(0.8, 1.0),(1.0, 1.2)]])
    topo = squaregrid(2, 3, 20.0, antiparallel=true)
    model, y, z, q = GndStr.make_master(inst, topo,
                                        SCIP.Optimizer(display_verblevel=0))

    JuMP.optimize!(model)
	status = JuMP.termination_status(model)
    @test status == :Optimal

    ysol = JuMP.value.(y)
    zsol = JuMP.value.(z)
    qsol = JuMP.value.(q)

    # shortest tree is obvious:
    @test ysol[11] ≈ 1.0
    @test ysol[13] ≈ 1.0
    @test ysol[4] ≈ 1.0
    @test sum(ysol) ≈ 3.0    # all others 0

    # only the smallest diameter is chosen
    @test zsol[11,1] ≈ 1.0
    @test zsol[13,1] ≈ 1.0
    @test zsol[4,1] ≈ 1.0
    @test sum(zsol) ≈ 3.0    # all others 0

    # uniq flow solution
    @test qsol[11] ≈ 50.0
    @test qsol[13] ≈ 20.0
    @test qsol[4] ≈ 30.0
    @test sum(qsol) ≈ 100.0  # all others 0
end

@testset "solve subproblem: ground structure, discrete diameters" begin
    #       7    9      even arc numbers for
    #   () - d2 - ()    reversed arcs
    #   /1   /3   /5
    #  s1 - () - d1
    #    11   13
    solver = GLPK.Optimizer()

    inst = Instance([Node(0,0), Node(40,0), Node(20, 20)],
                    [-50, 20, 30],
                    fill(Bounds(60,80), 3),
                    [Diameter(t...) for t in [(0.8, 1.0),(1.0, 1.2)]])
    topo = squaregrid(2, 3, 20.0, antiparallel=true)

    nnodes = 2*3
    narcs = 2*7
    ndiam = length(inst.diameters)

    # recreate solution from master problem
    zsol = fill(false, (narcs, ndiam))
    zsol[11,1] = true
    zsol[13,1] = true
    zsol[4,1] = true

    qsol = fill(0.0, narcs)
    qsol[11] = 50.0
    qsol[13] = 20.0
    qsol[4] = 30.0

    @testset "feasible subproblem" begin
        cand = GndStr.CandSol(zsol, qsol, fill(0.0, narcs))
        model, π, Δl, Δu, ploss, plb, pub = GndStr.make_sub(inst, topo, cand, solver)
        JuMP.optimize!(model)
	    status = JuMP.termination_status(model)
        @test status == :Optimal

        # primal solution
        πsol = JuMP.value.(π)
        Δlsol = JuMP.value.(Δl)
        Δusol = JuMP.value.(Δu)

        @test Δlsol ≈ zeros(nnodes)
        @test Δusol ≈ zeros(nnodes)

        # dual solution
        μsol = getdual(ploss)
        λlsol = getdual(plb)
        λusol = getdual(pub)

        @test μsol ≈ zeros(length(ploss))
        @test λlsol ≈ zeros(nnodes)
        @test λusol ≈ zeros(nnodes)
    end

    @testset "infeasible subproblem" begin
        cand = GndStr.CandSol(zsol, 10 * qsol, fill(0.0, narcs)) # scaled
        model, π, Δl, Δu, ploss, plb, pub = GndStr.make_sub(inst, topo, cand, solver)
        JuMP.optimize!(model)
	    status = JuMP.termination_status(model)
        @test status == :Optimal

        # primal solution
        πsol = JuMP.value.(π)
        Δlsol = JuMP.value.(Δl)
        Δusol = JuMP.value.(Δu)

        # complementarity
        @test Δlsol .* Δusol ≈ zeros(nnodes)

        # innodes should not have slack
        terms = [2, 3, 6]
        innodes = filter(v -> !(v in terms),  1:nnodes)
        @test Δlsol[innodes] ≈ zeros(length(innodes))
        @test Δusol[innodes] ≈ zeros(length(innodes))

        # at least two vertices have slack
        slack = Δlsol + Δusol
        @test sum(slack[terms] .> 0) >= 2

        # dual solution, TODO: fix sign of dual multipliers
        μsol = abs.(getdual(ploss))
        λlsol = abs.(getdual(plb))
        λusol = abs.(getdual(pub))

        # at least two bounds active
        @test sum(λlsol[terms] .> 0) >= 1
        @test sum(λusol[terms] .> 0) >= 1

        # at least one path active
        @test sum(μsol .> 0) >= 2

        @test λlsol[innodes] ≈ zeros(length(innodes))
        @test λusol[innodes] ≈ zeros(length(innodes))
    end
end

@testset "compare relaxation and exact for subproblem" begin
    #     __s1__  __s2
    #   t3      t4
    solver = GLPK.Optimizer()

    nodes = [Node(100,0), Node(300,0), Node(0,0), Node(200,0)]
    arcs = [Arc(1,3), Arc(1,4), Arc(2,4)]
    topo = Topology(nodes, arcs)

    demand = [-600, -400, 400, 600]
    bounds = fill(Bounds(40, 80), 4)
    diams = [Diameter(1.0, 1.0), Diameter(2.0, 2.0)]
    inst = Instance(nodes, demand, bounds, diams, ploss_coeff_nice)

    zsol = [true false; true false; true false]
    qsol = [400, 200, 400]
    cand = GndStr.CandSol(zsol, qsol, qsol.^2)

    @testset "solving the exact subproblem" begin
        model, π, Δl, Δu, ploss, plb, pub =
            GndStr.make_sub(inst, topo, cand, solver, relaxed=false)
        JuMP.optimize!(model)
	    status = JuMP.termination_status(model)
        @test status == :Optimal
        @test JuMP.objective_value(model) ≈ 3600
    end

    @testset "solving the relaxation" begin
        model, π, Δl, Δu, ploss, plb, pub =
            GndStr.make_sub(inst, topo, cand, solver, relaxed=true)
        JuMP.optimize!(model)
	    status = JuMP.termination_status(model)
        @test status == :Optimal
        @test JuMP.objective_value(model) ≈ 0
    end
end

@testset "run GBD iterations" begin
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

        result = optimize(inst, topo, GndStr.IterGBD(
            SCIP.Optimizer(display_verblevel=0),
            GLPK.Optimizer()))
        @test result.status == :Optimal

        zsol = result.solution.zsol
        @test zsol[4,1] == true
        @test zsol[11,1] == true
        @test zsol[13,1] == true
        @test sum(zsol) == 3

        @test result.dualbound ≈ 67.0
        @test result.niter == 1
    end

    @testset "medium flow: difficult instance" begin
        inst = Instance(nodes, 5*demand, bounds, diams)

        result = optimize(inst, topo, GndStr.IterGBD(
            SCIP.Optimizer(display_verblevel=0),
            GLPK.Optimizer()))
        @test result.status == :Optimal

        zsol = result.solution.zsol
        @test zsol[4,1] == true
        @test zsol[11,2] == true
        @test zsol[13,1] == true
        @test sum(zsol) == 3

        @test result.dualbound ≈ 72.0
        @test result.niter == 3
    end

    @testset "high flow: iteration limit instance" begin
        inst = Instance(nodes, 30*demand, bounds, diams)

        result = optimize(inst, topo, GndStr.IterGBD(
            SCIP.Optimizer(display_verblevel=0),
            GLPK.Optimizer(),
            maxiter=3))
        @test result.status == :UserLimit
        @test result.solution == nothing
        @test result.dualbound ≈ 156.0
        @test result.niter == 3
    end

    @testset "high flow: time limit instance" begin
        inst = Instance(nodes, 30*demand, bounds, diams)

        result = optimize(inst, topo, GndStr.IterGBD(
            SCIP.Optimizer(display_verblevel=0),
            GLPK.Optimizer(),
            timelimit=5.0))
        @test result.status == :UserLimit
        @test result.solution == nothing
    end

    @testset "high flow on triangle: infeasible" begin
        inst3 = Instance([Node(0,0), Node(50,0)],
                         20*[-50, 50],
                         [Bounds(60,80), Bounds(60,80)],
                         [Diameter(t...) for t in [(0.8, 1.0),(1.0, 1.2)]])
        topo3 = Topology([Node(0,0), Node(50,0), Node(30, 40)],
                         [Arc(1,3), Arc(1,2), Arc(2,3)])

        result = optimize(inst3, topo3, GndStr.IterGBD(
            SCIP.Optimizer(display_verblevel=0),
            GLPK.Optimizer()))
        @test result.status == :Infeasible
        @test result.solution == nothing
        @test result.dualbound == Inf
        @test result.niter == 2
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

        result = optimize(inst, topo, GndStr.IterGBD(
            SCIP.Optimizer(display_verblevel=0),
            GLPK.Optimizer()))
        @test result.status == :Optimal
        zsol = result.solution.zsol
        @test sum(zsol[:,2]) == 1
        @test sum(zsol[:,1]) == 2 # solution not quite uniqe
        @test result.dualbound ≈ 400.0
        @test result.niter == 2
    end

    @testset "difficult instance with disconnected candidates" begin
        inst = Instance([Node(0,0), Node(100,0), Node(200,100)],
                        [800, -900, 100],
                        fill(Bounds(40,80), 3),
                        [Diameter(t...) for t in [(1.0, 1.0), (2.0, 3.2)]],
                        ploss_coeff_nice)
        topo = squaregrid(2, 3, 100.0, antiparallel=true)

        # trigger the cuts for disconnected candidate
        solver = GndStr.IterGBD(SCIP.Optimizer(display_verblevel=0),
                                GLPK.Optimizer(),
                                addnogoods=true, addcritpath=false)
        result = optimize(inst, topo, solver)
        @test result.status == :Optimal

        zsol = result.solution.zsol
        @test sum(zsol) == 3

        @test result.dualbound ≈ 520.0
    end

end

@testset "Linear overestimation of supremum terms" begin
    values = Float64[i*j - 1 for i=1:4, j=1:4]
    values = 2*tril(values) - triu(values)

    cand_i, cand_j = 3, 2
    a, b, c = GndStr.linear_overest(values, cand_i, cand_j, GLPK.Optimizer())
    @test size(a) == (4,)
    @test size(b) == (4,)
    @test size(c) == ()

    @test a[cand_i] + b[cand_j] + c == values[cand_i, cand_j]
end

@testset "Solve semimaster without z vars" begin
    #       7    9      even arc numbers for
    #   () - d2 - ()    reversed arcs
    #   /1   /3   /5
    #  s1 - () - d1
    #    11   13
    inst = Instance([Node(0,0), Node(40,0), Node(20, 20)],
                    [-50, 20, 30],
                    fill(Bounds(60,80), 3),
                    [Diameter(t...) for t in [(0.8, 1.0),(1.0, 1.2)]])
    topo = squaregrid(2, 3, 20.0, antiparallel=true)
    model, y, q = GndStr.make_semimaster(inst, topo,
                                         SCIP.Optimizer(display_verblevel=0))

    JuMP.optimize!(model)
	status = JuMP.termination_status(model)
    @test status == :Optimal

    ysol = JuMP.value.(y)
    qsol = JuMP.value.(q)

    # shortest tree is obvious:
    @test ysol[11] ≈ 1.0
    @test ysol[13] ≈ 1.0
    @test ysol[4] ≈ 1.0
    @test sum(ysol) ≈ 3.0 # all others 0

    # uniq flow solution
    @test qsol[11] ≈ 50.0
    @test qsol[13] ≈ 20.0
    @test qsol[4] ≈ 30.0
    @test sum(qsol) ≈ 100.0 # all others 0
end

@testset "Solve semisubproblem with free z vars" begin
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

    zsol = fill(false, length(topo.arcs), length(diams))
    zsol[4,1]  = true
    zsol[11,1] = true
    zsol[13,1] = true
    qsol = fill(0.0, length(topo.arcs))
    qsol[[4, 11, 13]] = [30.0, 50.0, 20.0]

    @testset "low flow: very easy instance" begin
        factor = 1.0
        inst = Instance(nodes, factor*demand, bounds, diams)

        sol = GndStr.CandSol(zsol, factor*qsol, qsol.^2)
        model, candarcs, z = GndStr.make_semisub(
            inst, topo, sol, SCIP.Optimizer(display_verblevel=0))
        @test length(candarcs) == 3
        @test size(z) == (3, 2)

        JuMP.optimize!(model)
	    status = JuMP.termination_status(model)
        @test status == :Optimal

        znew = fill(false, length(topo.arcs), length(diams))
        znew[candarcs,:] = (JuMP.value.(z) .> 0.5)
        @test znew[4,1] == true
        @test znew[11,1] == true
        @test znew[13,1] == true
        @test sum(znew) == 3
    end

    @testset "medium flow: difficult instance" begin
        factor = 5.0
        inst = Instance(nodes, factor*demand, bounds, diams)

        sol = GndStr.CandSol(zsol, factor*qsol, qsol.^2)
        model, candarcs, z = GndStr.make_semisub(
            inst, topo, sol, SCIP.Optimizer(display_verblevel=0))
        @test length(candarcs) == 3
        @test size(z) == (3, 2)

        JuMP.optimize!(model)
	    status = JuMP.termination_status(model)
        @test status == :Optimal

        znew = fill(false, length(topo.arcs), length(diams))
        znew[candarcs,:] = (JuMP.value.(z) .> 0.5)
        @test znew[4,1] == true
        @test znew[11,2] == true
        @test znew[13,1] == true
        @test sum(znew) == 3

    end

    @testset "high flow on triangle: infeasible" begin
        inst = Instance([Node(0,0), Node(50,0)],
                        [-1000, 1000],
                        [Bounds(60,80), Bounds(60,80)],
                        [Diameter(t...) for t in [(0.8, 1.0),(1.0, 1.2)]])
        topo = Topology([Node(0,0), Node(50,0), Node(30, 40)],
                        [Arc(1,3), Arc(1,2), Arc(2,3)])

        zsol = fill(false, 3, 2)
        zsol[1,1] = true
        qsol = 1000*[1.0, 0.0, 0.0]
        sol = GndStr.CandSol(zsol, qsol, qsol.^2)
        model, candarcs, z = GndStr.make_semisub(
            inst, topo, sol, SCIP.Optimizer(display_verblevel=0))
        @test length(candarcs) == 1
        @test size(z) == (1, 2)

        JuMP.optimize!(model)
	    status = JuMP.termination_status(model)
        @test status == :Infeasible
    end

end

@testset "Solve semi decomposition with nogoods on y" begin
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
        inst = Instance(nodes, 1 * demand, bounds, diams)
        result = optimize(inst, topo, GndStr.IterTopo(
            SCIP.Optimizer(display_verblevel=0),
            SCIP.Optimizer(display_verblevel=0)))
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
        result = optimize(inst, topo, GndStr.IterTopo(
            SCIP.Optimizer(display_verblevel=0),
            SCIP.Optimizer(display_verblevel=0)))
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

        result = optimize(inst, topo, GndStr.IterTopo(
            SCIP.Optimizer(display_verblevel=0),
            SCIP.Optimizer(display_verblevel=0))
        @test result.status == :Infeasible
    end
end
