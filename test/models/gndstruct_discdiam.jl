import PipeLayout: squaregrid, CandSol, gndstruct_discdiam_master, gndstruct_discdiam_sub
using JuMP

facts("solve master problem (ground structure, discrete diameters)") do
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
    model, y, z, q = gndstruct_discdiam_master(inst, topo)

    status = solve(model)
    @fact status --> :Optimal

    ysol = getvalue(y)
    zsol = getvalue(z)
    qsol = getvalue(q)

    # shortest tree is obvious:
    @fact ysol[11] --> roughly(1.0)
    @fact ysol[13] --> roughly(1.0)
    @fact ysol[4] --> roughly(1.0)
    @fact sum(ysol) --> roughly(3.0) # all others 0

    # only the smallest diameter is chosen
    @fact zsol[11,1] --> roughly(1.0)
    @fact zsol[13,1] --> roughly(1.0)
    @fact zsol[4,1] --> roughly(1.0)
    @fact sum(zsol) --> roughly(3.0) # all others 0

    # uniq flow solution
    @fact qsol[11] --> roughly(50.0)
    @fact qsol[13] --> roughly(20.0)
    @fact qsol[4] --> roughly(30.0)
    @fact sum(qsol) --> roughly(100.0) # all others 0
end

facts("solve subproblem (ground structure, discrete diameters)") do
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

    context("feasible subproblem") do
        cand = CandSol(zsol, qsol)
        model, π, Δl, Δu = gndstruct_discdiam_sub(inst, topo, cand)
        status = solve(model)
        @fact status --> :Optimal

        πsol = getvalue(π)
        Δlsol = getvalue(Δl)
        Δusol = getvalue(Δu)

        @fact Δlsol --> roughly(zeros(nnodes))
        @fact Δusol --> roughly(zeros(nnodes))
    end

    context("infeasible subproblem") do
        cand = CandSol(zsol, 10 * qsol) # scaled
        model, π, Δl, Δu = gndstruct_discdiam_sub(inst, topo, cand)
        status = solve(model)
        @fact status --> :Optimal

        πsol = getvalue(π)
        Δlsol = getvalue(Δl)
        Δusol = getvalue(Δu)

        # complementarity
        @fact Δlsol .* Δusol --> roughly(zeros(nnodes))

        # innodes should not have slack
        terms = [2, 3, 6]
        isterm(v) = v in terms
        innodes = filter(not(isterm) , 1:nnodes)
        @fact Δlsol[innodes] --> roughly(zeros(3))
        @fact Δusol[innodes] --> roughly(zeros(3))

        # at least two vertices have slack
        slack = Δlsol + Δusol
        @fact sum(slack[terms] .> 0) --> greater_than_or_equal(2)
    end
end
