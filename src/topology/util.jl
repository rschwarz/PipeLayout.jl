using LightGraphs: Graph, DiGraph, add_edge!, neighbors

export arcindex, antiparallelindex, is_tree, find_cycle, pipelengths

"Build a LightGraphs.DiGraph from a Topology"
function digraph_from_topology(topology::Topology)
    n = length(topology.nodes)
    dg = DiGraph(n)
    for a in topology.arcs
        add_edge!(dg, a.tail, a.head)
    end
    dg
end

"Map to indices of arcs in topology."
function arcindex(topology::Topology)
    index = Dict{Arc,Int}()
    for (i,arc) in enumerate(topology.arcs)
        index[arc] = i
    end
    index
end

"""
Find indices of anti-parallel arcs.

Is 0 for arcs without anti-parallel arc.
"""
function antiparallelindex(topology::Topology)
    arcidx = arcindex(topology)
    result = fill(0, length(arcidx))
    for (a, (tail, head)) in enumerate(topology.arcs)
        anti = Arc(head, tail)
        if haskey(arcidx, anti)
            result[a] = arcidx[anti]
        end
    end
    result
end

"Checks whether given topology is a tree as undirected graph."
function is_tree(topology::Topology)
    n = length(topology.nodes)
    m = length(topology.arcs)
    m == n-1 && is_connected(digraph_from_topology(topology))
end

"Finds (undirected) cycle."
function find_cycle(topology::Topology)
    nnodes = length(topology.nodes)
    digraph = digraph_from_topology(topology)
    graph = Graph(digraph)

    # assume that graph is connected, so starting from single root is OK

    # initialize DFS
    root = 1
    stack = [(root, root)] # collect unvisited nodes with parent
    parent = fill(0, nnodes) # 0 means not visited

    # run search
    while length(stack) > 0 # not empty
        from, current = pop!(stack)
        parent[current] = from
        for to in neighbors(graph, current)
            if to == from
                continue # no antiparallel move
            elseif parent[to] == 0 # not visited
                push!(stack, (current, to))
            else # found cycle, run it in reverse
                path = [Arc(to, current)]
                backtrack = current
                while backtrack ≠ to
                    push!(path, Arc(backtrack, parent[backtrack]))
                    backtrack = parent[backtrack]
                end
                return path
            end
        end
    end

    # found no cycle
    return []
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
