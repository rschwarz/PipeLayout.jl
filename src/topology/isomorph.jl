using LightGraphs: Graph, neighbors, degree, nv, ne

export are_isomorphic

"find unlabeled nodes with at most one unlabeled neighbor"
function find_next_nodes_to_label(graph::Graph, labels::Vector{Int})
    n = nv(graph)
    @assert n == length(labels)
    nolabel = (1:n)[labels .== 0]
    filter(v -> sum(labels[neighbors(graph, v)] .== 0) ≤ 1, nolabel)
end

"combine neighbor labels into tentative labels"
function tentative_labels(graph::Graph, labels::Vector{Int}, nodes::Vector{Int})
    @assert all(labels[nodes] .== 0)
    tent = Tuple{Vector{Int},Int}[]
    for u in nodes
        nblabels = filter(l -> l ≠ 0, labels[neighbors(graph, u)])
        sort!(nblabels)
        push!(tent, (nblabels, u))
    end
    sort!(tent)
    tent
end

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

    # first label all the leaf nodes
    leaves1 = [n for n=1:nv(g1) if degree(g1, n) == 1]
    leaves2 = [n for n=1:nv(g2) if degree(g2, n) == 1]
    if length(leaves1) ≠ length(leaves2)
        return false
    end
    label1[leaves1] = label
    label2[leaves2] = label
    label += 1

    # add more labels until difference or isomorphism is found
    while true
        # find next nodes to give labels to
        next1 = find_next_nodes_to_label(g1, label1)
        next2 = find_next_nodes_to_label(g2, label2)
        if length(next1) ≠ length(next2)
            return false
        elseif length(next1) == 0
            return true
        end

        # try labels from multiset of neighbors' labels
        tent1 = tentative_labels(g1, label1, next1)
        tent2 = tentative_labels(g2, label2, next2)
        if any(t1[1] ≠ t2[1] for (t1, t2) in zip(tent1, tent2))
            return false
        end

        # actually assign the label to smallest of tentatives
        smallest = tent1[1][1]
        label1[[t[2] for t in tent1 if t[1] == smallest]] = label
        label2[[t[2] for t in tent2 if t[1] == smallest]] = label
        label += 1
    end
end

