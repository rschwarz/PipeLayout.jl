using PipeLayout
import PipeLayout: has_geosteiner, euclidean_steiner_tree

facts("compute ESMT using GeoSteiner") do
    # skip these tests if GeoSteiner is not installed
    has_geosteiner() || return

    nodes = [Node(0,0), Node(2,0), Node(1,1)]
    topo = euclidean_steiner_tree(nodes)
    @fact length(topo.nodes) --> 4
    @fact topo.nodes[1] --> nodes[1]
    @fact topo.nodes[2] --> nodes[2]
    @fact topo.nodes[3] --> nodes[3]
    @fact topo.nodes[4].x --> roughly(1.0)
    @fact topo.nodes[4].y --> roughly(0.57735, atol=1e-5)

    @fact length(topo.arcs) --> 3
    for arc in topo.arcs
        @fact findfirst(arc, 4) --> anyof(1,2)
    end
end
