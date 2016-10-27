using PipeLayout
using FactCheck

include("helpers.jl")

include("gasphysics.jl")
include("util.jl")
include("pwl.jl")

include("topology/geosteiner.jl")
include("topology/grid.jl")
include("topology/mst.jl")
include("topology/util.jl")

include("flow.jl")

include("models/pipedimensioning.jl")
include("models/gndstruct_discdiam.jl")

if Pkg.installed("SCIP") != nothing
    include("models/junctionlocation.jl")
end

FactCheck.exitstatus()
