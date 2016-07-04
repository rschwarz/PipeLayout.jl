using DataStructures: DisjointSets, union!, in_same_set, num_groups

"Compute array of edges (Nodes in columns), sorted by length"
function sorted_edges(nodes::Vector{Node})
    edges = hcat(combinations(nodes, 2)...)
    sortcols(edges, by=edge -> begin
             u, v = edge
             (u.x - v.x)^2 + (u.y - v.y)^2
             end)
end

"""
Compute a minimum spanning tree using Kruskal's algorithm.

Returns a vector of arcs with integer indices for nodes.
"""
function min_spanning_tree(nodes::Vector{Node})
    # MAYBE, use a Delaunay triangulation instead?
    n = length(nodes)

    edges = sorted_edges(nodes)
    @assert size(edges, 1) == 2 # column-wise

    nodeidx = [nodes[i] => i for i in 1:n]

    tree_edges = Arc[]
    sizehint!(tree_edges, n-1)

    sets = DisjointSets{Node}(nodes)
    for e in 1:size(edges, 2)
        u, v = edges[:, e]
        if !in_same_set(sets, u, v)
            union!(sets, u, v)
            push!(tree_edges, Arc(nodeidx[u], nodeidx[v]))
        end
        num_groups(sets) == 1 && break # spanning tree
    end

    tree_edges
end

"Build a Topology from an Instance's nodes, using MST"
function topology_from_mst(nodes::Vector{Node})
    Topology(nodes, min_spanning_tree(nodes))
end
topology_from_mst(instance::Instance) = topology_from_mst(instance.nodes)
