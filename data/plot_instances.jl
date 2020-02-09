using JSON
using PipeLayout

open(ARGS[1]) do f
    for line in eachline(f)
        key = strip(line)
        println(key)
        inst = PipeLayout.read_instance(".", key)
        topo = PipeLayout.read_topology(".", key)

        scene = PipeLayout.empty_scene()
        draw!(scene, topo)
        # TODO: draw node demand
        save("$key.png", scene)
    end
end
