import PipeLayout: dist

@testset "digraph and incidence" begin
    # FST example:
    #  2      4
    #   \    /
    #    5--6
    #   /    \
    #  1      3
    nodes = [Node(i,i) for i in 1:6]
    arcs = [Arc(5, 6), Arc(5,1), Arc(2, 5), Arc(3, 6), Arc(6, 4)]
    topo = Topology(nodes, arcs)

    N = incidence(topo)
    @test size(N) == (6, 5)
    @test N == [ 0.0  1.0  0.0  0.0  0.0
                 0.0  0.0 -1.0  0.0  0.0
                 0.0  0.0  0.0 -1.0  0.0
                 0.0  0.0  0.0  0.0  1.0
                -1.0 -1.0  1.0  0.0  0.0
                 1.0  0.0  0.0  1.0 -1.0]
end

@testset "check tree topology" begin
    nodes = [Node(i,i) for i in 1:3]
    arcs = [Arc(1,2), Arc(1, 3), Arc(2, 3)]

    @testset "disconnected forest" begin
        @test !is_tree(Topology(nodes, arcs[1:1]))
    end

    @testset "path" begin
        @test is_tree(Topology(nodes, arcs[1:2]))
    end

    @testset "triangle" begin
        @test !is_tree(Topology(nodes, arcs[1:3]))
    end
end

@testset "check cycle finding" begin
    nodes = [Node(i,i^2) for i in 1:4]
    arcs = [Arc(1,2), Arc(2,3), Arc(3,4), Arc(2,4)]

    @test find_cycle(Topology(nodes, arcs[1:3])) == []

    result = find_cycle(Topology(nodes, arcs))
    @test length(result) == 3
    @test result[1].head == result[2].tail
    @test result[2].head == result[3].tail
    @test result[3].head == result[1].tail
end

@testset "check node distance and pipe lengths" begin
    u, v, w = Node(0,0), Node(0,3), Node(4,0)

    @test dist(u, u) ≈ 0
    @test dist(u, v) ≈ 3
    @test dist(u, w) ≈ 4
    @test dist(v, w) ≈ 5

    arcs = [Arc(1,2), Arc(1, 3), Arc(2, 3)]
    topo = Topology([u, v, w], arcs)

    @test pipelengths(topo) ≈ [3, 4, 5]
end

@testset "check arcindex and antiparallelindex" begin
    nodes = [Node(0,0), Node(0,1), Node(1,0)]
    arcs = [Arc(1,2), Arc(1,3), Arc(2,1)]
    topo = Topology(nodes, arcs)

    n, m = 3, 3

    arcidx = arcindex(topo)
    @test length(arcidx) == m
    @test arcidx[Arc(1,2)] == 1
    @test arcidx[Arc(2,1)] == 3
    @test haskey(arcidx, Arc(1,3)) == true
    @test haskey(arcidx, Arc(3,1)) == false

    antiidx = antiparallelindex(topo)
    @test length(antiidx) == m
    @test antiidx[1] == 3
    @test antiidx[2] == 0
    @test antiidx[3] == 1
end
