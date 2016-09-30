using PipeLayout

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
    for nterminals in [3, 5, 7, 11]
        key = @sprintf "%s-t%02d"  gskey nterminals
        terminals[key] = select_subset(topo.nodes, nterminals)
    end
end

# for each set of terminals (and related topo), create demand situations
#  [ ] random distribution (source vs sink)
#  [ ] same, but single source
#  [ ] hard distribution: maximize transport momentum
#  [ ] for each normalized demand: several flow scalings

# uniform pressure bounds
#  [ ] one for each boundary node

# diameter data
#  [ ] full data and smaller subset

# create instances by combination
#  [ ] come up with uniq naming
#  [ ] write instances to files (json)
