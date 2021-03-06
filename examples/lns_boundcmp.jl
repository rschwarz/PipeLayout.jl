# we use warmstart to hint at a solution that has better objective value in pipe
# dimensioning than the dual bound of junction location?!

using MathOptInterface
using JuMP
using SCIP

using PipeLayout
using PipeLayout.JuncLoc
using PipeLayout.PipeDim
import PipeLayout: ploss_coeff_nice

const MOI = MathOptInterface

# instance data
pipelength = 100.0
width = height = pipelength / sqrt(2)
terminals = [Node(0,0), Node(width,height), Node(2width,0), Node(3width,height)]
demand = [400, -600, 600, -400]
bounds = fill(Bounds(40, 80), 4)
diams = map(Diameter, [1, 2], [1, 2])
inst = Instance(terminals, demand, bounds, diams, ploss_coeff_nice)

_scip = JuMP.optimizer_with_attributes(SCIP.Optimizer, "limits/gap" => 0.001)

function pipedim()
    println("### Pipe Dimensioning LP model, path topology ###")

    # topology for path
    path = Topology(terminals, [Arc(2,1), Arc(2,3), Arc(4,3)])

    # instead of calling `optimize` we do the steps manually
    solver = PipeDim.LP(_scip)
    model, π, l = PipeDim.make_model(inst, path, solver.lpsolver)
    JuMP.optimize!(model)

    println(model)
    @show JuMP.value.(l)
    @show JuMP.value.(π)
end

function juncloc()
    println("### Junction Location NLP model, FST with fixed locations ###")

    # topology data for FST
    junctions = [Node(0,-1), Node(0,-2)]
    nodes = vcat(terminals, junctions)
    fstw = Topology(nodes, [Arc(5,6), Arc(5,1), Arc(2,5), Arc(3,6), Arc(6,4)])

    # instead of calling `optimize` we do the steps manually
    solver = JuncLoc.NLP(_scip)
    model, x, y, L, l, π = JuncLoc.make_nlp(inst, fstw, solver)

    # fix values for the Steiner node positions
    # 1st steiner is at 2nd term
    @constraint(model, x[5] == inst.nodes[2].x)
    @constraint(model, y[5] == inst.nodes[2].y)
    # 2nd steiner is at 3rd term
    @constraint(model, x[6] == inst.nodes[3].x)
    @constraint(model, y[6] == inst.nodes[3].y)

    JuMP.optimize!(model)

    println(model)
    @show JuMP.value.(L)
    @show JuMP.value.(l)
    @show JuMP.value.(π)
end

pipedim()
println()
juncloc()

# demands = [400, -600, 600, -400]
#
# LP:                flow   flow^2
# 1: 1 --> 2  (a)    400     40000
# 2: 2 --> 3  (b)    200    160000
# 3: 3 --> 4  (c)    400     40000
#
#      2   4
#     / \ /
#    1   3
#
#
# NLP:                flow   flow squared
# 1: 5 --> 6 (b)      200       40000
# 2: 5 --> 1 (a)      400      160000
# 3: 2 --> 5 (-)      600      360000
# 4: 3 --> 6 (-)     -600      360000
# 5: 6 --> 4 (c)     -400      160000
#
#       2     _4
#     _ 5 -- 6
#    1       3
#
