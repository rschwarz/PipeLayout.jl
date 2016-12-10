import PipeLayout: sorted_edges, min_spanning_tree, topology_from_mst, digraph_from_topology
using PipeLayout
using LightGraphs

@testset "compute MST from nodes" begin
    coords = [0 1 0; 0 0 2]
    nodes = Node[Node(coords[:,j]...) for j in 1:size(coords, 2)]
    const n = length(nodes)
    @test n == 3

    @testset "test sorted edges" begin
        edges = sorted_edges(nodes)
        m = n*(n-1)/2
        @test m == 3
        @test size(edges, 1) == 2
        @test size(edges, 2) == m
        @test edges[:, 1] == Node[nodes[1], nodes[2]]
        @test edges[:, 2] == Node[nodes[1], nodes[3]]
        @test edges[:, 3] == Node[nodes[2], nodes[3]]
    end

    @testset "look at resulting MST" begin
        const n = length(nodes)
        const m = n - 1 # tree

        tree = nodes |> topology_from_mst |> digraph_from_topology
        @test isa(tree, DiGraph)
        @test is_directed(tree)
        @test nv(tree) == n
        @test ne(tree) == m
        @test is_connected(tree)
        @test has_edge(tree, 1, 2)
        @test has_edge(tree, 1, 3)
    end
end
