using Plots: plot, scatter!, text, with
using Colors: colormap

immutable Net
    n::Int64
    m::Int64
    posx::Vector{Float64}
    posy::Vector{Float64}
    edges::Array{Int64, 2}

    function Net(posx::Vector{Float64}, posy::Vector{Float64}, edges::Array{Int64, 2})
        n = length(posx)
        m = size(edges, 2)

        @assert length(posx) == length(posy)
        @assert size(edges, 1) == 2
        @assert all(k -> 1 <= k <= n, edges)

        new(n, m, posx, posy, edges)
    end
end

"""Compute an interval with a relative margin around all values"""
function bbox(values; margin=0.05)
    lb, ub = extrema(values)
    range = lb == ub ? 1.0 : ub - lb
    lb - margin*range, ub + margin*range
end


"""Low level drawing method for given node positions and styling.

Vertices are assumed to by identified by 1:n and edges is a 2×m array
referencing the vertices. The style can be given as value (for all edges),
or as array (per edge).
"""
function draw(posx, posy, edges;
              edgecolor=:gray, edgewidth=2,
              nodecolor=:orange, nodesize=12)
    @assert size(posx) == size(posy)
    @assert size(posx, 2) == 1
    @assert size(edges, 1) == 2
    n = size(posx, 1)
    m = size(edges, 2)

    with(leg=false, grid=false, aspect_ratio=1, xlim=bbox(posx), ylim=bbox(posx)) do
        # draw edges
        plot(posx[edges], posy[edges], color=edgecolor, linewidth=edgewidth)

        # draw nodes
        scatter!(posx, posy, marker=(nodesize, nodecolor),
                 seriesann=map(x -> text(x, nodesize-3), 1:n))
    end
end

"scale and shift numbers to interval [lb=0, ub=1]"
function normalize{T<:Real}(data::Array{T}; lb=zero(T), ub=one(T))
    @assert ub > lb
    mi, ma = extrema(data)
    range = ma == mi ? 1.0 : ma - mi
    (data - mi) / range * (ub - lb) + lb
end

"map numbers to colors"
function map2color(data::Array{Float64}; ncolors::Int=20, cmap="Blues")
    data = normalize(data, lb=0.0, ub=1.0 - 1E-9)

    colors = colormap(cmap, ncolors)
    indices = round(Int, ncolors*data + 1, RoundDown)
    colors[indices]
end

"drawing of networks, with mapping of attributes to color and size"
function draw(net::Net;
              edgecolor=:gray, edgecmap="Blues",
              edgewidth=2, edgebds=[2, 4],
              nodecolor=:orange, nodecmap="Reds",
              nodesize=12)
    edgecolor = (typeof(edgecolor) == Array{Float64, 2}
                 ? map2color(edgecolor, cmap=edgecmap)
                 : edgecolor)
    edgewidth = (typeof(edgewidth) == Array{Float64, 2}
                 ? normalize(edgewidth, lb=edgebds[1], ub=edgebds[2])
                 : edgewidth)
    nodecolor = (typeof(nodecolor) == Array{Float64, 2}
                 ? map2color(nodecolor, cmap=nodecmap)
                 : nodecolor)

    draw(net.posx, net.posy, net.edges,
         edgecolor=edgecolor, edgewidth=edgewidth,
         nodecolor=nodecolor, nodesize=nodesize)
end

function draw(topo::Topology)
    n, m = length(topo.nodes), length(topo.arcs)
    nodepos = reshape(reinterpret(Float64, topo.nodes), (2,n))
    posx, posy = nodepos[1,:]', nodepos[2,:]'
    arcs = reshape(reinterpret(Int, topo.arcs), (2,m))
    draw(posx, posy, arcs)
end
