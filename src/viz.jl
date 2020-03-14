using AbstractPlotting: Scene, scatter!, linesegments!

function empty_scene(; kwargs...)
    return Scene(; kwargs..., show_axis=false, scale_plot=false)
end

function draw!(scene, topo::Topology; markersize=2)
    # arcs (bottom layer)
    positions = Array(topo.nodes[vcat(collect.(topo.arcs)...)])
    linesegments!(scene, positions)

    # nodes (on top)
    positions = topo.nodes
    scatter!(scene, positions, markersize=markersize)
end
