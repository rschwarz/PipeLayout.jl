# TODO: move this graph code somewhere else!
using LightGraphs: DiGraph, add_edge!
using DataStructures: DisjointSets, union!, in_same_set, num_groups

immutable PipeDimInstance
    instance::Instance
    tree::DiGraph
end

"Compute array of edges (Nodes in columns), sorted by length"
function sorted_edges(nodes::Vector{Node})
    edges = hcat(combinations(nodes, 2)...)
    sortcols(edges, by=edge -> begin
             u, v = edge
             (u.x - v.x)^2 + (u.y - v.y)^2
             end)
end


"""Compute a minimum spanning tree using Kruskal's algorithm

Returns a LightGraphs.DiGraphs with integer indices for nodes.
"""
function min_spanning_tree(nodes::Vector{Node})
    # TODO, use a Delaunay triangulation instead?
    n = length(nodes)
    tree = DiGraph(n)

    edges = sorted_edges(nodes)
    @assert size(edges, 1) == 2 # column-wise

    nodeidx = [nodes[i] => i for i in 1:n]

    sets = DisjointSets{Node}(nodes)
    for e in 1:size(edges, 2)
        u, v = edges[:, e]
        if !in_same_set(sets, u, v)
            union!(sets, u, v)
            add_edge!(tree, nodeidx[u], nodeidx[v])
        end
        num_groups(sets) == 1 && break # spanning tree
    end

    tree
end

"Build a PipeDimInstance from an Instance, using MST for the topology"
function instance_from_mst(instance::Instance)
    PipeDimInstance(instance, min_spanning_tree(instance))
end

# TODO: actually build the model, think about naming
function build_pd_model()
    m = Model()
end
