using AbstractPlotting: Scene, scatter!, linesegments!
using CairoMakie

function empty_scene(; kwargs...)
    return Scene(; kwargs..., show_axis=false, scale_plot=false)
end

function draw!(scene, topo::Topology)
    # arcs (bottom layer)
    positions = Array(topo.nodes[vcat(topo.arcs...)])
    linesegments!(scene, positions)

    # nodes (on top)
    positions = topo.nodes
    scatter!(scene, positions)
end
