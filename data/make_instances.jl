using PipeLayout

# create topologies for ground structures
#  [x] square meshes of increasing sizes
#  [ ] random (or sobolev) points with voronoi triangulation
#  [ ] grids with holes
ground_structures = Topology[]
append!(ground_structures,
        [squaregrid(args..., antiparallel=true) for args in
         [(6, 7, 15.0), (9, 9, 10.0), (13, 13, 7.5), (19, 19, 5.0)]])
# TODO: also create names for these, and store them in Dict?

# for each topology: select subsets of nodes for terminals
#  [ ] randomly
#  [ ] greedy farthest neighbor

# for each set of terminals (and related topo), create demand situations
#  [ ] random distribution (source vs sink)
#  [ ] hard distribution: maximize transport momentum
#  [ ] for each normalized demand: several flow scalings

# uniform pressure bounds
#  [ ] one for each boundary node

# diameter data
#  [ ] full data and smaller subset

# create instances by combination
#  [ ] come up with uniq naming
#  [ ] write instances to files (json)
