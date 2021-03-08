activenv("diss")
using Revise

using Random

using PipeLayout

using AbstractPlotting
using CairoMakie
using FileIO


function data(N=7)
    Random.seed!(0);
    WIDTH = 500
    HEIGHT = 300

    x = round.(0.95 * WIDTH * rand(N))
    y = round.(0.95 * HEIGHT * rand(N))
    points = collect([x y]')
end

function draw(path, trimesh)
    scene = PipeLayout.empty_scene()
    topo = convert(Topology, trimesh)
    PipeLayout.draw!(scene, topo)
    save(File(format"PNG", path), scene)
    nothing
end

points = data()
del = PipeLayout.delaunay_triangulation(points)
draw("del.png", del)

six = PipeLayout.refine_sixths(del)
draw("six.png", six)

six2 = PipeLayout.refine_sixths(del, minimum_angle=15)
draw("six2.png", six2)
s = PipeLayout.TriangleSwitches(minimum_angle=20, maximum_area=2000)
ref = PipeLayout.refine(six2, s)
draw("ref.png", ref)
