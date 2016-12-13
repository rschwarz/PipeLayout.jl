#
# test data
#
"create a small instance (and topology) with a single pipe"
function single_pipe(;length=100.0, flow=200.0)
    nodes = [Node(0,0), Node(length,0)]
    demand = [-flow, flow]
    press = fill(Bounds(60.0, 80.0), size(nodes))
    diams = [Diameter(t...) for t in [(0.8, 1.0),(1.0, 1.2),(1.2, 1.5)]]
    inst = Instance(nodes, demand, press, diams)
    topo = Topology(nodes, [Arc(1,2)])
    inst, topo
end
