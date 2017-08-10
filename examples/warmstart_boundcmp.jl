# we use warmstart to hint at a solution that has better objective value in pipe
# dimensioning than the dual bound of junction location?!

using JuMP
using SCIP

using PipeLayout
using PipeLayout.JuncLoc
import PipeLayout: ploss_coeff_nice

# instance data
pipelength = 100.0
width = height = pipelength / sqrt(2)
terminals = [Node(0,0), Node(width,height), Node(2width,0), Node(3width,height)]
demand = [400, -600, 600, -400]
bounds = fill(Bounds(40, 80), 4)
diams = map(Diameter, [1, 2], [1, 2])
inst = Instance(terminals, demand, bounds, diams, ploss_coeff_nice)

# topology data for FST
junctions = [Node(0,-1), Node(0,-2)]
nodes = vcat(terminals, junctions)
fstw = Topology(nodes, [Arc(5,6), Arc(5,1), Arc(2,5), Arc(3,6), Arc(6,4)])

# instead of calling `optimize` we do the steps manually
solver = JuncLoc.NLP(SCIPSolver("limits/gap", 0.01,
                                "heuristics/completesol/maxunkownrate", 1.0))
model, x, y, L, l, π = PipeLayout.JuncLoc.make_nlp(inst, fstw, solver)

# add dummy binary to support incomplete warmstart
@variable(model, dummybin, Bin)

# set start values for the Steiner node positions
# 1st steiner is at 2nd term
setvalue(x[5], inst.nodes[2].x)
setvalue(y[5], inst.nodes[2].y)
# 2nd steiner is at 3rd term
setvalue(x[6], inst.nodes[3].x)
setvalue(y[6], inst.nodes[3].y)

status = solve(model)
