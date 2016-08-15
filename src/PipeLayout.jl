module PipeLayout

include("instances.jl")
# --> deserialization.jl
# --> random.jl

include("util.jl")
include("gasphysics.jl")

# graph methods and topologies
include("topology/util.jl")
include("topology/geosteiner.jl")
include("topology/grid.jl")
include("topology/mst.jl")

include("flow.jl") # needs topology

# optimization models and approaches
module PipeDimensioning
importall ..PipeLayout
include("models/pipedimensioning.jl")
end

include("models/gndstruct_discdiam.jl")

# visualization
include("draw.jl")

end
