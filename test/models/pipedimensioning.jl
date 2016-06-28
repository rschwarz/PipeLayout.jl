import PipeLayout: sorted_edges

facts("compute sorted edges from nodes") do
    coords = [0 1 0;
              0 0 2]
    nodes = [Node(coords[:,j]...) for j in 1:size(coords, 2)]
    n = length(nodes)
    @fact n --> 3

    edges = sorted_edges(nodes)
    m = n*(n-1)/2
    @fact m --> 3
    @fact size(edges, 1) --> 2
    @fact size(edges, 2) --> m
    @fact edges[:, 1] --> Node[nodes[1], nodes[2]]
    @fact edges[:, 2] --> Node[nodes[1], nodes[3]]
    @fact edges[:, 3] --> Node[nodes[2], nodes[3]]
end
