using Printf
using Random

using JSON
using PipeLayout

# set random seed for reproducible results
Random.seed!(23);

# create positions for terminals randomly
const WIDTH = 120 - 1
const HEIGHT = 90 - 1
randnode() = Node(rand(0:WIDTH), rand(0:HEIGHT))
terminals = Dict{String, Vector{Node}}()
for nterminals in [3, 5, 7, 9, 11]
    for variant in ("a", "b")
        key = @sprintf "t%02d%s" nterminals variant
        terminals[key] = [randnode() for _ in 1:nterminals]
    end
end

# create ground structures
# TODO

# for each set of terminals, create demand situations
#  [x] random distribution (source vs sink)
#  [x] same, but single source
#  [X] for each normalized demand: several flow scalings
demands = Dict{String, Vector{Float64}}()
for (tkey, nodes) in terminals
    nnodes = length(nodes)
    nodeidx = collect(1:nnodes)

    # collect assignments to source/sink role
    sourcesinks = []

    # single source, at first position
    push!(sourcesinks, (nodeidx[1:1], nodeidx[2:end]))

    # half and half (less sources)
    if nnodes > 3
        half = div(nnodes, 2)
        push!(sourcesinks, (nodeidx[1:half], nodeidx[half+1:end]))
    end

    for (sources, sinks) in sourcesinks
        for variant in ["a", "b"]
            dist = PipeLayout.randdemand(nnodes, sources, sinks, 1.0)
            for flow in [25, 50, 100, 200]
                key = @sprintf "%s.%02d_%02d_%s_%03d" tkey length(sources) length(sinks) variant flow
                demands[key] = flow * dist
            end
        end
    end
end

# uniform pressure bounds
#  [x] one for each boundary node
pressure_bounds = Bounds(40.0, 80.0)

# diameter data
#  [x] full data and smaller subset
_diams = [0.4,  0.6,  0.8, 1.0,  1.2]
_costs = [0.68, 0.91, 1.2, 1.55, 1.96]
diameters = Dict{String, Vector{Diameter}}()
diameters["full"] = [Diameter(d,c) for (d,c) in zip(_diams, _costs)]
diameters["some"] = diameters["full"][1:2:end]

# create instances by combination
#  [x] come up with uniq naming
#  [x] write instances to files (json)
instances = String[]
for (key, demand) in demands
    gs, trm, dmd = split(key, ".")
    topology = ground_structures[gs]
    nodes = terminals["$gs.$trm"]
    pressure = [pressure_bounds for n in nodes]

    for (diamkey, diams) in diameters
        instkey = "$key.$diamkey"
        instance = Instance(nodes, demand, pressure, diams, ploss_coeff_nice)
        push!(instances, instkey)

        open("$instkey.instance.json", "w") do f
            JSON.print(f, instance, 2)
        end
        open("$instkey.topology.json", "w") do f
            JSON.print(f, topology, 2)
        end
    end
end

sort!(instances)
open("instances.txt", "w") do f
    write(f, join(instances, "\n"))
end
