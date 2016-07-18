using LightGraphs: DiGraph, add_edge!

function digraph_from_topology(topology::Topology)
    n = length(topology.nodes)
    dg = DiGraph(n)
    for a in topology.arcs
        add_edge!(dg, a.tail, a.head)
    end
    dg
end

function arcindex(topology::Topology)
    index = Dict{Arc,Int}()
    for (i,arc) in enumerate(topology.arcs)
        index[arc] = i
    end
    index
end

function is_tree(topology::Topology)
    n = length(topology.nodes)
    m = length(topology.arcs)
    m == n-1 && is_connected(digraph_from_topology(topology))
end

"Euclidean distance of two nodes"
function dist(u::Node, v::Node)
    norm(u - v)
end

"Compute lengths of all pipes"
function pipelengths(topo::Topology)
    lengths = Float64[]
    sizehint!(lengths, length(topo.arcs))
    for arc in topo.arcs
        tail = topo.nodes[arc.tail]
        head = topo.nodes[arc.head]
        push!(lengths, dist(tail, head))
    end
    lengths
end
