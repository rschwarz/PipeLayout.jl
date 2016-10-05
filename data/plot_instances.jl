using JSON
using PipeLayout
using Plots

pyplot()

open(ARGS[1]) do f
    for line in eachline(f)
        key = strip(line)
        println(key)
        inst = PipeLayout.read_instance(".", key)
        topo = PipeLayout.read_topology(".", key)

        plot(inst)
        plot!(topo)
        savefig("$key.png")
    end
end
