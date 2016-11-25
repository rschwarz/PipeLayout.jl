facts("check tree isomorphy") do
    context("fails on non-trees") do
        nodes = [Node(i,i) for i in 1:3]
        arcs = [Arc(1,2), Arc(1, 3), Arc(2, 3)]

        triangle = Topology(nodes, arcs)
        @fact_throws are_isomorphic(triangle, triangle)

        forest = Topology(nodes, arcs[1:1])
        @fact_throws are_isomorphic(forest, forest)
    end

    context("test on trees") do
        nodes = [Node(i,i) for i in 1:5]
        path1 = Topology(nodes, [Arc(i,i+1) for i=1:4])
        path2 = Topology(nodes, [Arc(1,3), Arc(3,5), Arc(5,2), Arc(2,4)])
        star1 = Topology(nodes, [Arc(1,i) for i=2:5])
        star2 = Topology(nodes, [Arc(5,i) for i=1:4])

        @fact are_isomorphic(path1, path1) --> true
        @fact are_isomorphic(path1, path2) --> true
        @fact are_isomorphic(star1, star1) --> true
        @fact are_isomorphic(star1, star2) --> true

        @fact are_isomorphic(path1, star1) --> false
        @fact are_isomorphic(star1, path2) --> false
    end
end

facts("test steiner tree enumeration") do
    nclass3, repr3 = enumerate_steiner(3)
    @fact length(nclass3) --> 1
    @fact nclass3[1] --> 1
    @fact length(repr3) --> sum(nclass3)

    nclass4, repr4 = enumerate_steiner(4)
    @fact length(nclass4) --> 2
    @fact nclass4 --> [1,1]
    @fact length(repr4) --> sum(nclass4)

    nclass5, repr5 = enumerate_steiner(5)
    @fact length(nclass5) --> 3
    @fact nclass5 --> [1,1,1]
    @fact length(repr5) --> sum(nclass5)

    nclass6, repr6 = enumerate_steiner(6)
    @fact length(nclass6) --> 4
    @fact nclass6 --> [1,1,1,2]
    @fact length(repr6) --> sum(nclass6)
end

facts("test FST labeling") do
    nodes = [Node(i,i) for i in 1:6]

    # single steiner node
    Y = Topology(nodes[1:4], [Arc(1,4), Arc(2,4), Arc(3,4)])
    @fact label_fst(Y) --> (2, 3)

    # two steiner nodes
    H = Topology(nodes, [Arc(1,5), Arc(2,5), Arc(3,6),Arc(4,6), Arc(5,6)])
    @fact label_fst(H) --> (2, (3, 4))
end
