import PipeLayout: is_tree

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
