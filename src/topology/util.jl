using LightGraphs: DiGraph, add_edge!

function digraph_from_topology(topology::Topology)
    n = length(topology.nodes)
    dg = DiGraph(n)
    for a in topology.arcs
        add_edge!(dg, a.tail, a.head)
    end
    dg
end

function is_tree(topology::Topology)
    n = length(topology.nodes)
    m = length(topology.arcs)
    m == n-1 && is_connected(digraph_from_topology(topology))
end
