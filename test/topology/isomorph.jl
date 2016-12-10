@testset "check tree isomorphy" begin
    @testset "fails on non-trees" begin
        nodes = [Node(i,i) for i in 1:3]
        arcs = [Arc(1,2), Arc(1, 3), Arc(2, 3)]

        triangle = Topology(nodes, arcs)
        @test_throws AssertionError are_isomorphic(triangle, triangle)

        forest = Topology(nodes, arcs[1:1])
        @test_throws AssertionError are_isomorphic(forest, forest)
    end

    @testset "test on trees" begin
        nodes = [Node(i,i) for i in 1:5]
        path1 = Topology(nodes, [Arc(i,i+1) for i=1:4])
        path2 = Topology(nodes, [Arc(1,3), Arc(3,5), Arc(5,2), Arc(2,4)])
        star1 = Topology(nodes, [Arc(1,i) for i=2:5])
        star2 = Topology(nodes, [Arc(5,i) for i=1:4])

        @test are_isomorphic(path1, path1) == true
        @test are_isomorphic(path1, path2) == true
        @test are_isomorphic(star1, star1) == true
        @test are_isomorphic(star1, star2) == true

        @test are_isomorphic(path1, star1) == false
        @test are_isomorphic(star1, path2) == false
    end
end

@testset "test steiner tree enumeration" begin
    nclass3, repr3 = enumerate_steiner(3)
    @test length(nclass3) == 1
    @test nclass3[1] == 1
    @test length(repr3) == sum(nclass3)

    nclass4, repr4 = enumerate_steiner(4)
    @test length(nclass4) == 2
    @test nclass4 == [1,1]
    @test length(repr4) == sum(nclass4)

    nclass5, repr5 = enumerate_steiner(5)
    @test length(nclass5) == 3
    @test nclass5 == [1,1,1]
    @test length(repr5) == sum(nclass5)

    nclass6, repr6 = enumerate_steiner(6)
    @test length(nclass6) == 4
    @test nclass6 == [1,1,1,2]
    @test length(repr6) == sum(nclass6)
end

@testset "test FST labeling" begin
    nodes = [Node(i,i) for i in 1:6]

    # single steiner node
    Y = Topology(nodes[1:4], [Arc(1,4), Arc(2,4), Arc(3,4)])
    @test label_fst(Y) == (2, 3)

    # two steiner nodes
    H = Topology(nodes, [Arc(1,5), Arc(2,5), Arc(3,6),Arc(4,6), Arc(5,6)])
    @test label_fst(H) == (2, (3, 4))
end

@testset "test FST enumeration" begin
    terminals = [Node(2,3), Node(4,5), Node(6,7), Node(8,9)]

    @test_throws AssertionError enumerate_fst([])
    @test_throws AssertionError enumerate_fst(terminals[1:1])
    @test_throws AssertionError enumerate_fst(terminals[1:2])

    trees3 = enumerate_fst(terminals[1:3])
    @test length(trees3) == 1

    trees4 = enumerate_fst(terminals[1:4])
    @test length(trees4) == 1 + 2
end
