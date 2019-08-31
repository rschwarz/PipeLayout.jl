using PipeLayout
using PipeLayout.GndStr

using JuMP
using GLPK

nodes = [Node(0,0), Node(45,0), Node(25, 22)]
demand = 10*[-50, 20, 30]
bounds = fill(Bounds(60,80), 3)
diams = [Diameter(t...) for t in [(0.8, 1.0),(1.0, 1.2)]]

topo = Topology([Node(t...) for t in [(0,22), (0,0), (25,22), (25,0), (45,22), (45,0)]],
                [Arc(t...) for t in [(1,2), (2,1), (3,4), (4,3), (5,6), (6,5),
                                     (1,3), (3,1), (3,5), (5,3),
                                     (2,4), (4,2), (4,6), (6,4)]])

inst = Instance(nodes, demand, bounds, diams, ploss_coeff_nice)

_glpk = JuMP.with_optimizer(GLPK.Optimizer)
solver = GndStr.IterGBD(_glpk, _glpk,
                        debug=true, addnogoods=false, maxiter=10)
result = optimize(inst, topo, solver);

@show result.status
@show result.dualbound
@show result.niter
