import PipeLayout:  dist
using PipeLayout

facts("build grid topology with 4 neighbors") do
    vert, hori, degree, width = 3, 4, 4, 5.0
    grid = RegGrid{vert,hori,degree}(width)
    topo = gridtopology(grid)
    nodes = topo.nodes
    arcs = topo.arcs

    narcs = 2*((vert - 1)*hori + vert*(hori - 1))
    @fact length(nodes) --> vert*hori
    @fact length(arcs) --> narcs
    for arc in arcs
        tail, head = arc
        @fact dist(nodes[tail], nodes[head]) --> roughly(width)
    end

    # bottom left node should be at origin
    nodemat = reshape(nodes, (vert, hori))
    @fact nodemat[vert, 1].x --> roughly(0.0)
    @fact nodemat[vert, 1].y --> roughly(0.0)
end

facts("build grid topology with 8 neighbors") do
    vert, hori, degree, width = 3, 4, 8, 5.0
    grid = RegGrid{vert,hori,degree}(width)
    topo = gridtopology(grid)
    nodes = topo.nodes
    arcs = topo.arcs

    narcs = 2*((vert - 1)*hori + vert*(hori - 1) + 2*(vert - 1)*(hori - 1))
    @fact length(nodes) --> vert*hori
    @fact length(arcs) --> narcs

    diag = sqrt(2.0) * width
    for arc in arcs
        tail, head = arc
        @fact dist(nodes[tail], nodes[head]) --> anyof(roughly(width), roughly(diag))
    end

    # bottom left node should be at origin
    nodemat = reshape(nodes, (vert, hori))
    @fact nodemat[vert, 1].x --> roughly(0.0)
    @fact nodemat[vert, 1].y --> roughly(0.0)
end
