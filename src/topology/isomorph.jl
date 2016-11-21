using LightGraphs: Graph, neighbors, degree

export are_isomorphic

"""
Are two given (undirected) trees isomorphic?

Implemented by uncritical copy from Fampa et al.: A specialized branch-and-bound
algorithm for the Euclidean Steiner tree problem in n-space", algorithm 2.
"""
function are_isomorphic(tree1::Topology, tree2::Topology)
    @assert is_tree(tree1) && is_tree(tree2)
    g1 = Graph(digraph_from_topology(tree1))
    g2 = Graph(digraph_from_topology(tree2))

    if nv(g1) != nv(g2) || ne(g1) != ne(g2)
        return false
    end

    label1, label2 = fill(0, nv(g1)), fill(0, nv(g2))
    label = 1

    leaves1 = [n for n=1:nv(g1) if degree(g1, n) == 1]
    leaves2 = [n for n=1:nv(g2) if degree(g2, n) == 1]

    # TODO...
end
