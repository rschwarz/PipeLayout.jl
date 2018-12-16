module PipeLayout

include("instances.jl")
# --> deserialization.jl
# --> random.jl
include("types.jl")

include("util.jl")
include("gasphysics.jl")
include("pwl.jl")
include("diameter.jl")

# graph methods and topologies
include("topology/util.jl")
include("topology/geosteiner.jl")
include("topology/grid.jl")
include("topology/mst.jl")
include("topology/isomorph.jl")

include("flow.jl") # needs topology

# optimization models and approaches
module PipeDim
using ..PipeLayout
include("models/pipedimensioning.jl")
end

module JuncLoc
using ..PipeLayout
include("models/junctionlocation/main.jl")
end

module GndStr
using ..PipeLayout
include("models/gndstruct_discdiam/main.jl")
end

# visualization via Plots recipes
include("plot_recipes.jl")

end
