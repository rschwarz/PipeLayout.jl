import PipeLayout: uniq_flow, flow_path_decomp, reorient_fwdflow

@testset "compute unique flow on trees" begin
    nodes = [Node(i,i) for i in 1:6]
    n = length(nodes)

    @testset "everything fine (H net)" begin
        arcs = [Arc(1,3), Arc(2,3), Arc(3,4), Arc(4,5), Arc(4,6)]
        d = [-2.0, -4.0, 0.0, 0.0, 3.0, 3.0]
        q = uniq_flow(Topology(nodes, arcs), d)
        @test length(q) == length(arcs)
        @test q ≈ [-d[1], -d[2], d[5] + d[6], d[5], d[6]]
    end

    @testset "everything fine (H net with Steiner nodes)" begin
        arcs = [Arc(1,3), Arc(2,3), Arc(3,4), Arc(4,5), Arc(4,6)]
        topo = Topology(nodes, arcs)

        d = [-2.0, -4.0, 3.0, 3.0]
        terms = nodes[[1, 2, 5, 6]]
        inst = Instance(terms, d, fill(Bounds(1, 2), 4), [Diameter(1, 1)])

        q = uniq_flow(inst, topo)
        @test length(q) == length(arcs)
        @test q ≈ [-d[1], -d[2], d[3] + d[4], d[3], d[4]]
    end

    @testset "everything fine (FST with 4 terms + 2 stein)" begin
        arcs = [Arc(5,6), Arc(5,1), Arc(2,5), Arc(3,6), Arc(6,4)]
        topo = Topology(nodes, arcs)
        d = [400.0, -600.0, 600.0, -400.0, 0.0, 0.0]

        q = uniq_flow(topo, d)
        @test length(q) == length(arcs)
        @test q ≈ [200, -400, 600, -600, 400]
    end

    @testset "wrong topology" begin
        arcs = [Arc(1,2), Arc(2,3), Arc(3,4), Arc(4,5), Arc(5,6), Arc(1,6)]
        d = [-2.0, -4.0, 0.0, 0.0, 3.0, 3.0]
        @test_throws ArgumentError uniq_flow(Topology(nodes, arcs), d)
    end

    @testset "unbalanced demand" begin
        arcs = [Arc(1,3), Arc(2,3), Arc(3,4), Arc(4,5), Arc(4,6)]
        d = fill(1.0, 6)
        @test_throws ArgumentError uniq_flow(Topology(nodes, arcs), d)
    end

    @testset "dimension mismatch" begin
        arcs = [Arc(1,2), Arc(2,3), Arc(3,4), Arc(4,5), Arc(5,6)]
        d = [-2.0, -4.0, 0.0, 0.0, 3.0, 3.0, 1.0]
        @test_throws ArgumentError uniq_flow(Topology(nodes, arcs), d)
    end
end

@testset "reorient arcs for forward flow" begin
    # star with arcs going to center
    nodes = [Node(i,i) for i=1:4]
    arcs = [Arc(i,4) for i=1:3]
    topo = Topology(nodes, arcs)

    @testset "nothing changes for zero flow" begin
        rotopo = reorient_fwdflow(topo, fill(0.0, 4))
        @test rotopo.nodes == topo.nodes
        @test rotopo.arcs == topo.arcs
    end

    @testset "change arcs leading to sink" begin
        rotopo = reorient_fwdflow(topo, [-1.0, 1.0, 0.0, 0.0])
        @test rotopo.nodes == topo.nodes
        @test rotopo.arcs[1] == topo.arcs[1] # source
        @test rotopo.arcs[2].tail == topo.arcs[2].head # sink
        @test rotopo.arcs[2].head == topo.arcs[2].tail # sink
        @test rotopo.arcs[3] == topo.arcs[3] # innode
    end
end

@testset "decompose positive arc flow in paths" begin
    @testset "on a simple tree" begin
        topo = Topology([Node(0,0), Node(1,0), Node(2,0), Node(1,1)],
                        [Arc(1,2), Arc(2,3), Arc(2,4)])
        arcflow = [3.0, 2.0, 1.0]

        paths, pathflows = flow_path_decomp(topo, arcflow)
        @test length(paths) == 2
        @test length(pathflows) == 2
        @test paths[1] == [Arc(1,2), Arc(2,3)]
        @test pathflows[1] ≈ 2.0
        @test paths[2] == [Arc(1,2), Arc(2,4)]
        @test pathflows[2] ≈ 1.0
    end

    @testset "on a grid with flow on subtree" begin
        # ()-1-()-6-
        # 2    7
        # ()-4-()-10-
        grid = RegGrid{2,3,4}(50.0)
        topo = gridtopology(grid)
        nnodes, narcs = length(topo.nodes), length(topo.arcs)
        arcflow = fill(0.0, narcs)
        arcflow[1] = 3.0
        arcflow[6] = 2.0
        arcflow[7] = 1.0

        paths, pathflows = flow_path_decomp(topo, arcflow)
        @test length(paths) == 2
        @test length(pathflows) == 2
        @test paths[1] == [topo.arcs[1], topo.arcs[6]]
        @test pathflows[1] ≈ 2.0
        @test paths[2] == [topo.arcs[1], topo.arcs[7]]
        @test pathflows[2] ≈ 1.0
    end
end
