import PipeLayout:  dist
using PipeLayout

@testset "build grid topology with 4 neighbors" begin
    vert, hori, degree, width = 3, 4, 4, 5.0
    grid = RegGrid{vert,hori,degree}(width)
    topo = gridtopology(grid)
    nodes = topo.nodes
    arcs = topo.arcs

    narcs = 2*((vert - 1)*hori + vert*(hori - 1))
    @test length(nodes) == vert*hori
    @test length(arcs) == narcs
    for arc in arcs
        tail, head = arc
        @test dist(nodes[tail], nodes[head]) ≈ width
    end

    # bottom left node should be at origin
    nodemat = reshape(nodes, (vert, hori))
    @test nodemat[vert, 1].x ≈ 0.0
    @test nodemat[vert, 1].y ≈ 0.0
end

@testset "build grid topology with 8 neighbors" begin
    vert, hori, degree, width = 3, 4, 8, 5.0
    grid = RegGrid{vert,hori,degree}(width)
    topo = gridtopology(grid)
    nodes = topo.nodes
    arcs = topo.arcs

    narcs = 2*((vert - 1)*hori + vert*(hori - 1) + 2*(vert - 1)*(hori - 1))
    @test length(nodes) == vert*hori
    @test length(arcs) == narcs

    diag = sqrt(2.0) * width
    for arc in arcs
        tail, head = arc
        d = dist(nodes[tail], nodes[head])
        @test  d ≈ width || d ≈ diag
    end

    # bottom left node should be at origin
    nodemat = reshape(nodes, (vert, hori))
    @test nodemat[vert, 1].x ≈ 0.0
    @test nodemat[vert, 1].y ≈ 0.0
end
