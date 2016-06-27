# TODO: move this graph code somewhere else!
using LightGraphs

immutable PipeDimInstance
    instance::Instance
    tree::DiGraph
end

"Compute array of edges (as columns), sorted by length"
function sorted_edges(nodes::Vector{Node})
    edges = hcat(combinations(nodes, 2)...)
    sortcols(edges, by=edge -> begin
             u, v = edge
             (u.x - v.x)^2 + (u.y - v.y)^2
             end)
end


"Compute a minimum spanning tree using Kruskal's algorithm"
function min_spanning_tree(nodes::Vector{Node})
    # TODO, use a Delaunay triangulation instead?
    edges = sorted_edges(nodes)
    # now, use DisjointSets from DataStructures
    nothing
end

# TODO: actually build the model, think about naming
function build_pd_model()
    m = Model()
end
