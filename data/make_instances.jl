using PipeLayout
import PipeLayout: randdemand

# set random seed for reproducible results
srand(23)

# create topologies for ground structures
#  [x] square meshes of increasing sizes
#  [ ] other meshes (with crossing diagonals? cf. Jakob's thesis)
#  [ ] random (or sobolev) points with voronoi triangulation
#  [ ] grids with holes
ground_structures = Dict{String, Topology}()
for args in [(6, 7, 15.0), (9, 9, 10.0), (13, 13, 7.5), (19, 19, 5.0)]
    key = @sprintf "m%02d_n%02d_l%04.1f" args...
    ground_structures[key] = squaregrid(args..., antiparallel=true)
end

# for each topology: select subsets of nodes for terminals
#  [x] randomly
#  [ ] greedy farthest neighbor
terminals = Dict{String, Vector{Node}}()
for (gskey, topo) in ground_structures
    for nterminals in [3, 5, 7]
        for variant in ["a", "b"]
            key = @sprintf "%s-t%02d_%s" gskey nterminals variant
            terminals[key] = select_subset(topo.nodes, nterminals)
        end
    end
end

# for each set of terminals (and related topo), create demand situations
#  [ ] random distribution (source vs sink)
#  [ ] same, but single source
#  [ ] hard distribution: maximize transport momentum
#  [ ] for each normalized demand: several flow scalings
demands = Dict{String, Vector{Float64}}()
for (tkey, nodes) in terminals
    nnodes = length(nodes)
    nodeidx = collect(1:nnodes)

    # single source, at first position
    sources, sinks = nodeidx[1:1], nodeidx[2:end]
    for variant in ["a", "b"]
        dist = randdemand(nnodes, sources, sinks, 1.0)
        for flow in [25, 50, 100, 200]
            key = @sprintf "%s-%02d_%02d_%s_%03d" tkey 1 length(sinks) variant flow
            demands[key] = flow * dist
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
diams_full = [Diameter(d,c) for (d,c) in zip(_diams, _costs)]
diams_some = diams_full[1:2:end]

# create instances by combination
#  [x] come up with uniq naming
#  [ ] write instances to files (json)
