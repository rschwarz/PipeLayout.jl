using Colors: colormap, RGB
using RecipesBase

"scale and shift numbers to interval [lb=0, ub=1]"
function normalize{T<:Real}(data::Array{T}; lb=zero(T), ub=one(T))
    @assert ub > lb
    mi, ma = extrema(data)
    range = ma == mi ? 1.0 : ma - mi
    (data - mi) / range * (ub - lb) + lb
end

"map numbers to colors"
function map2color(data::Array{Float64}; ncolors::Int=20, cmap="Blues")
    data = normalize(data, lb=0.1, ub=0.9 - 1E-9)

    colors = colormap(cmap, ncolors)
    indices = round(Int, ncolors*data + 1, RoundDown)
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
    nodes, Arc[]
end

@recipe function f(topo::PipeLayout.Topology)
    topo.nodes, topo.arcs
end
