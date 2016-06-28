import PipeLayout: sorted_edges, min_spanning_tree
using PipeLayout
using LightGraphs

facts("compute MST from nodes") do
    coords = [0 1 0;
              0 0 2]
    nodes = [Node(coords[:,j]...) for j in 1:size(coords, 2)]
    n = length(nodes)
    @fact n --> 3

    context("test sorted edges") do
        edges = sorted_edges(nodes)
        m = n*(n-1)/2
        @fact m --> 3
        @fact size(edges, 1) --> 2
        @fact size(edges, 2) --> m
        @fact edges[:, 1] --> Node[nodes[1], nodes[2]]
        @fact edges[:, 2] --> Node[nodes[1], nodes[3]]
        @fact edges[:, 3] --> Node[nodes[2], nodes[3]]
    end

    context("look at resulting MST") do
        tree = min_spanning_tree(nodes)
        @fact tree --> is_directed
        @fact nv(tree) --> 3
        @fact ne(tree) --> 2
        @fact tree --> is_connected
        @fact has_edge(tree, 1, 2) --> true
        @fact has_edge(tree, 1, 3) --> true
    end
end
