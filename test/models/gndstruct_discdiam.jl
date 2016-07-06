import PipeLayout: squaregrid, gndstruct_discdiam_master
using JuMP

facts("solve master problem") do
    #   () - d2 - ()
    #   /    /    /
    #  s1 - () - d1
    inst = Instance([Node(0,0), Node(40,0), Node(20, 20)],
                    [-50, 20, 30],
                    fill(Bounds(60,80), 3),
                    [Diameter(t...) for t in [(0.8, 1.0),(1.0, 1.2)]])
    topo = squaregrid(2, 3, 20.0, antiparallel=true)
    model, y, z, q = gndstruct_discdiam_master(inst, topo)

    # writeLP(model, "/home/bzfschwa/foo.lp")

    status = solve(model)
    @fact status --> :Optimal
end
