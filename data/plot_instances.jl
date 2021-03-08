using AbstractPlotting
using CairoMakie
using PipeLayout

function draw_instance(key)
    inst = PipeLayout.read_instance(".", key)
    topo = PipeLayout.read_topology(".", key)

    scene = PipeLayout.empty_scene()
    PipeLayout.draw!(scene, topo)
    # TODO: draw node demand
    save("$key.png", scene)
end

function draw_instances(filename)
    open(filename) do f
        for line in eachline(f)
            key = strip(line)
            draw_instance(key)
        end
    end
end

draw_instances(ARGS[1])
