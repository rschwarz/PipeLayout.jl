using Test
using PipeLayout
using Aqua

include("helpers.jl")

@testset "basics" begin
    include("gasphysics.jl")
    include("util.jl")
    include("pwl.jl")
    include("diameter.jl")
end

@testset "topology" begin
    include("topology/geosteiner.jl")
    include("topology/grid.jl")
    include("topology/mst.jl")
    include("topology/util.jl")
    include("topology/isomorph.jl")
end

@testset "flow" begin
    include("flow.jl")
end

@testset "pipe dimensioning" begin
    include("models/pipedimensioning.jl")
end

@testset "ground structure & discrete diameter" begin
    include("models/gndstruct_discdiam.jl")
end

@testset "junction location" begin
    include("models/junctionlocation.jl")
    include("models/junctionlocation_nlp.jl")
end

Aqua.test_all(PipeLayout)
