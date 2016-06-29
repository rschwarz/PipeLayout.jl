using LightGraphs: DiGraph, add_edge!

function digraph_from_topology(topology::Topology)
    n = length(topology.nodes)
    dg = DiGraph(n)
    for a in topology.arcs
        add_edge!(dg, a.tail, a.head)
    end
    dg
end
