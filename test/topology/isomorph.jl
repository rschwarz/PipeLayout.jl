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
