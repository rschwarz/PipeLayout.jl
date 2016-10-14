using JSON
using PipeLayout
import PipeLayout: ploss_coeff_nice, randdemand

# set random seed for reproducible results
srand(23)

# create topologies for ground structures
#  [x] square meshes of increasing sizes
#  [x] other meshes (with crossing diagonals? cf. Jakob's thesis)
#  [ ] random (or sobolev) points with voronoi triangulation
#  [ ] grids with holes
ground_structures = Dict{String, Topology}()
for args in [(6, 7, 18), (9, 9, 12), (13, 13, 9), (19, 19, 6)]
    m, n, w = args
    for deg in [4, 8]
        key = @sprintf "m%02d_n%02d_d%02d_l%02d" n m deg w
        ground_structures[key] = gridtopology(RegGrid{m,n,deg}(w))
    end
end

# for each topology: select subsets of nodes for terminals
#  [x] randomly
#  [ ] greedy farthest neighbor
terminals = Dict{String, Vector{Node}}()
for (gskey, topo) in ground_structures
    for nterminals in [3, 5, 7]
        key = @sprintf "%s.t%02d" gskey nterminals
        terminals[key] = select_subset(topo.nodes, nterminals)
    end
end

# for each set of terminals (and related topo), create demand situations
#  [x] random distribution (source vs sink)
#  [x] same, but single source
#  [ ] hard distribution: maximize transport momentum
#  [ ] for each normalized demand: several flow scalings
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
            dist = randdemand(nnodes, sources, sinks, 1.0)
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
