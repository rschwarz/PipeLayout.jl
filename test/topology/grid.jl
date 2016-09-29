import PipeLayout:  dist
using PipeLayout

facts("build grid topology") do
    vert, hori, width = 3, 4, 5.0
    topo = squaregrid(vert, hori, width)
    nodes = topo.nodes
    arcs = topo.arcs

    @fact length(nodes) --> vert*hori
    @fact length(arcs) --> (vert - 1)*hori + vert*(hori - 1)
    for arc in arcs
        tail, head = arc
        @fact dist(nodes[tail], nodes[head]) --> roughly(width)
    end

    # bottom left node should be at origin
    nodemat = reshape(nodes, (vert, hori))
    @fact nodemat[vert, 1].x --> roughly(0.0)
    @fact nodemat[vert, 1].y --> roughly(0.0)
end

facts("build grid topology with antiparallel arcs") do
    vert, hori, width = 3, 4, 5.0
    topo = squaregrid(vert, hori, width, antiparallel=true)
    nodes = topo.nodes
    arcs = topo.arcs

    narcs = 2*((vert - 1)*hori + vert*(hori - 1))
    @fact length(nodes) --> vert*hori
    @fact length(arcs) --> narcs
    for arc in arcs
        tail, head = arc
        @fact dist(nodes[tail], nodes[head]) --> roughly(width)
    end

    # anti-parallel arcs
    for a in 1:2:narcs
        fwarc, bwarc = arcs[a], arcs[a+1]
        @fact fwarc.tail --> bwarc.head
        @fact fwarc.head --> bwarc.tail
    end

    # bottom left node should be at origin
    nodemat = reshape(nodes, (vert, hori))
    @fact nodemat[vert, 1].x --> roughly(0.0)
    @fact nodemat[vert, 1].y --> roughly(0.0)
end
