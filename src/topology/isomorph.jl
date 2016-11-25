import Base:isless
using LightGraphs: Graph, neighbors, degree, nv, ne

export are_isomorphic, enumerate_steiner, label_fst

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

Reference: algorithm 2 from Fampa et al.: A specialized branch-and-bound
algorithm for the Euclidean Steiner tree problem in n-space.
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

"""
Enumerate all representative trees.

This considers only the trees interconnecting the Steiner nodes for full Steiner
trees of n terminals.

Following algorithm 1 in Fampa et al.: A specialized branch-and-bound algorithm
for the Euclidean Steiner tree problem in n-space.
"""
function enumerate_steiner(n::Int)
    nodes = [Node(i,i) for i=1:n-2]
    nclass = fill(0, n - 2) # by number of steiner nodes
    repr = Dict()

    # start with single node
    nclass[1] = 1
    repr[(1, 1)] = Topology(nodes[1:1], Arc[])

    # add more steiner nodes
    for j in 2:(n - 2)
        # look at trees of one size smaller
        for k in 1:nclass[j-1]
            t_base = repr[(j-1, k)]
            g_base = Graph(digraph_from_topology(t_base))

            # and all potential neighbors
            for i in 1:(j-1)
                # maximum degree
                if degree(g_base, i) == 3
                    continue
                end

                # connect new node here then for new tree
                t_new = Topology(nodes[1:j], vcat(t_base.arcs, [Arc(i,j)]))

                # avoid isomorphies
                if any(are_isomorphic(t_new, repr[(j,l)]) for l in 1:nclass[j])
                    continue
                end

                # add to collection
                nclass[j] += 1
                repr[(j, nclass[j])] = t_new
            end
        end
    end
    nclass, repr
end

# needed below
isless{T<:Int}(x::T, y::Tuple{T,T}) = true
isless{T<:Int}(x::Tuple{T,T}, y::T) = false

"""
Computes a labeling of a Full Steiner Tree.

This can be used to check for isomorphic trees (having same label).

Following algorithm 4 in Fampa et al.: A specialized branch-and-bound algorithm
for the Euclidean Steiner tree problem in n-space.
"""
function label_fst(fst::Topology)
    n = length(fst.nodes)
    g = Graph(digraph_from_topology(fst))
    terms = [v for v in 1:n if degree(g, v) == 1]
    stein = [v for v in 1:n if degree(g, v) == 3]
    @assert length(terms) + length(stein) == n

    label = Any[-1 for i = 1:n]
    label[terms] = terms

    # rooted binary tree with lexicographical ordering
    root = minimum(terms)
    current_level, next_level = [root], []
    bfs_edges = []
    visited = fill(false, n)
    # BFS
    while length(current_level) > 0
        visited[current_level] = true
        empty!(next_level)
        for v in current_level
            for w in neighbors(g, v)
                if !visited[w]
                    push!(next_level, w)
                    push!(bfs_edges, (v,w))
                end
            end
        end
        current_level = next_level[:] # copy!
    end

    # work up from the leaves
    reverse!(bfs_edges)
    for depth in 1:2:length(bfs_edges)-1
        t1, h1 = bfs_edges[depth]
        t2, h2 = bfs_edges[depth + 1]
        @assert t1 == t2
        @assert label[t1] == -1
        l1, l2 = label[h1], label[h2]
        label[t1] = l1 < l2 ? (l1, l2) : (l2, l1)
    end

    # all labeled
    @assert all(label .!= -1)

    # label of root's adjacent Steiner node is label for tree
    t, h = bfs_edges[end]
    @assert t == root
    label[h]
end
