import PipeLayout: squaregrid, gndstruct_discdiam_master
using JuMP

facts("solve master problem") do
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
