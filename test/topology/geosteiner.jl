using PipeLayout
import PipeLayout: has_geosteiner, euclidean_steiner_tree

if has_geosteiner()
    @testset "compute ESMT using GeoSteiner" begin
        nodes = [Node(0,0), Node(2,0), Node(1,1)]
        topo = euclidean_steiner_tree(nodes)
        @test length(topo.nodes) == 4
        @test topo.nodes[1] == nodes[1]
        @test topo.nodes[2] == nodes[2]
        @test topo.nodes[3] == nodes[3]
        @test topo.nodes[4].x ≈ 1.0
        @test topo.nodes[4].y ≈ 0.57735 atol=1e-5

        @test length(topo.arcs) == 3
        for arc in topo.arcs
            @test findfirst(arc, 4) in [1,2]
        end
    end
end
