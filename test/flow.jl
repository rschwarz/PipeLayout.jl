import PipeLayout: uniq_flow, flow_path_decomp

facts("compute unique flow on trees") do
    nodes = [Node(i,i) for i in 1:6]
    n = length(nodes)

    context("everything fine (H net)") do
        arcs = [Arc(1,3), Arc(2,3), Arc(3,4), Arc(4,5), Arc(4,6)]
        m = length(arcs)
        d = [-2.0, -4.0, 0.0, 0.0, 3.0, 3.0]
        q = uniq_flow(Topology(nodes, arcs), d)
        @fact length(q) --> m
        @fact q --> roughly([-d[1], -d[2], d[5] + d[6], d[5], d[6]])
    end

    context("wrong topology") do
        arcs = [Arc(1,2), Arc(2,3), Arc(3,4), Arc(4,5), Arc(5,6), Arc(1,6)]
        d = [-2.0, -4.0, 0.0, 0.0, 3.0, 3.0]
        @fact_throws ArgumentError uniq_flow(Topology(nodes, arcs), d)
    end

    context("unbalanced demand") do
        arcs = [Arc(1,3), Arc(2,3), Arc(3,4), Arc(4,5), Arc(4,6)]
        d = fill(1.0, 6)
        @fact_throws ArgumentError uniq_flow(Topology(nodes, arcs), d)
    end

    context("dimension mismatch") do
        arcs = [Arc(1,2), Arc(2,3), Arc(3,4), Arc(4,5), Arc(5,6)]
        d = [-2.0, -4.0, 0.0, 0.0, 3.0, 3.0, 1.0]
        @fact_throws ArgumentError uniq_flow(Topology(nodes, arcs), d)
    end
end

facts("decompose positive arc flow in paths") do
    context("on a simple tree") do
        topo = Topology([Node(0,0), Node(1,0), Node(2,0), Node(1,1)],
                        [Arc(1,2), Arc(2,3), Arc(2,4)])
        arcflow = [3.0, 2.0, 1.0]

        paths, pathflows = flow_path_decomp(topo, arcflow)
        @fact length(paths) --> 2
        @fact length(pathflows) --> 2
        @fact paths[1] --> [Arc(1,2), Arc(2,3)]
        @fact pathflows[1] --> roughly(2.0)
        @fact paths[2] --> [Arc(1,2), Arc(2,4)]
        @fact pathflows[2] --> roughly(1.0)
    end

    context("on a grid with flow on subtree") do
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
        @fact length(paths) --> 2
        @fact length(pathflows) --> 2
        @fact paths[1] --> [topo.arcs[1], topo.arcs[6]]
        @fact pathflows[1] --> roughly(2.0)
        @fact paths[2] --> [topo.arcs[1], topo.arcs[7]]
        @fact pathflows[2] --> roughly(1.0)
    end
end
