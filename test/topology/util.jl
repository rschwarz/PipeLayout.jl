import PipeLayout: is_tree, dist, pipelengths

facts("check tree topology") do
    nodes = [Node(i,i) for i in 1:3]
    arcs = [Arc(1,2), Arc(1, 3), Arc(2, 3)]

    context("disconnected forest") do
        @fact Topology(nodes, arcs[1:1]) --> not(is_tree)
    end

    context("path") do
        @fact Topology(nodes, arcs[1:2]) --> is_tree
    end

    context("triangle") do
        @fact Topology(nodes, arcs[1:3]) --> not(is_tree)
    end
end

facts("check node distance and pipe lengths") do
    u, v, w = Node(0,0), Node(0,3), Node(4,0)

    @fact dist(u, u) --> roughly(0)
    @fact dist(u, v) --> roughly(3)
    @fact dist(u, w) --> roughly(4)
    @fact dist(v, w) --> roughly(5)

    arcs = [Arc(1,2), Arc(1, 3), Arc(2, 3)]
    topo = Topology([u, v, w], arcs)

    @fact pipelengths(topo) --> roughly([3, 4, 5])
end
