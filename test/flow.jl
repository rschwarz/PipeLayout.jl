import PipeLayout: uniq_flow

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
