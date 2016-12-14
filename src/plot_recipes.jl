using Colors: colormap
using RecipesBase

"shift and scale numbers to interval [lb=0, ub=1]"
function shiftnscale{T<:Real}(data::Array{T}; lb=zero(T), ub=one(T))
    @assert ub > lb
    mi, ma = extrema(data)
    range = ma == mi ? 1.0 : ma - mi
    (data - mi) / range * (ub - lb) + lb
end

"map numbers to colors"
function map2color{T<:Real}(data::Array{T}; ncolors::Int=20, cmap="Blues")
    data = shiftnscale(data)

    # take two extra colors, one to leave out at each extremity, and another one
    # because the normalized data is in [0, 1] inclusive
    colors = colormap(cmap, ncolors + 3)
    indices = round(Int, ncolors*data + 2, RoundDown)
    colors[indices]
end

# a recipe to display topologies with arcs and nodes
@recipe function f(nodes::Array{PipeLayout.Node}, arcs::Array{PipeLayout.Arc})
    # prepare node data
    nodepos = reshape(reinterpret(Float64, nodes), (2, length(nodes)))'
    x = nodepos[:, 1]
    y = nodepos[:, 2]

    # prepare arc data
    arcidx = reshape(reinterpret(Int64, arcs), (2, length(arcs)))

    # first draw the arcs
    @series begin
        color --> :black
        x[arcidx], y[arcidx]
    end

    # then the nodes
    @series begin
        seriestype := :scatter
        aspect_ratio --> :equal
        markercolor --> :white
        markersize --> 3
        grid --> false
        foreground_color_border	--> :white
        legend --> false
        ticks --> nothing
        x, y
    end
end

# a recipe to display only nodes
@recipe function f(nodes::Array{PipeLayout.Node})
    markeralpha --> 0.5
    markercolor --> :red
    markersize --> 6
    nodes, Arc[]
end

# a recipe for topologies: just nodes and arcs
@recipe function f(topo::PipeLayout.Topology)
    topo.nodes, topo.arcs
end

# a recipe to display an instance: nodes with demand
@recipe function f(instance::PipeLayout.Instance)
    nodes = instance.nodes
    demand = instance.demand

    colors = fill(:gray, length(nodes))
    colors[demand .> 0] = :red
    colors[demand .< 0] = :blue

    sizes = sqrt(abs(demand))
    sizes = 10 * sizes / maximum(sizes) + 4

    markercolor --> colors
    markersize --> sizes
    nodes, Arc[]
end

# a recipe to display a candidate solution: active arcs with diameter choices
@recipe function f(topo::PipeLayout.Topology,
                   cand::PipeLayout.GndStr.CandSol)
    narcs, ndiams = size(cand.zsol)
    arcidx, diamidx = findn(cand.zsol)
    arcs = topo.arcs[arcidx]

    linewidth --> 3 * diamidx'
    linecolor --> map2color(diamidx', ncolors=ndiams, cmap="Greens")

    # need to redraw nodes, because arcs reference them
    markersize --> 0

    topo.nodes, arcs
end

# a recipe to display a pipe dimensioning solution: diameter choices
@recipe function f(inst::PipeLayout.Instance,
                   topo::PipeLayout.Topology,
                   sol::PipeLayout.PipeDim.Solution)
    narcs, ndiams = size(sol.lsol)
    equiv = [serialmerge(inst.diameters, sol.lsol[a,:]) for a=1:narcs]

    linewidth --> 8 * equiv'
    linecolor --> map2color(equiv', ncolors=ndiams, cmap="Greens")

    # need to redraw nodes, because arcs reference them
    markersize --> 0

    topo.nodes, topo.arcs
end

# a recipe to display a junction location solution: nodes & diameter choices
@recipe function f(inst::PipeLayout.Instance,
                   topo::PipeLayout.Topology,
                   sol::PipeLayout.JuncLoc.Solution)
    nodes = sol.nodes
    narcs, ndiams = size(sol.lsol)
    equiv = [serialmerge(inst.diameters, sol.lsol[a,:]) for a=1:narcs]

    linewidth --> 8 * equiv'
    linecolor --> map2color(equiv', ncolors=ndiams, cmap="Greens")

    # need to redraw nodes, because arcs reference them.
    # also: show junction nodes (not in Instance)
    markersize --> 8

    sol.nodes, topo.arcs
end
