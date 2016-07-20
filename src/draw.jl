using Plots: plot, plot!, scatter, scatter!, text, with
using Colors: colormap, RGB

export draw

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
function boundbox(values; margin=0.05)
    lb, ub = extrema(values)
    range = lb == ub ? 1.0 : ub - lb
    lb - margin*range, ub + margin*range
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
    data = normalize(data, lb=0.1, ub=0.9 - 1E-9)

    colors = colormap(cmap, ncolors)
    indices = round(Int, ncolors*data + 1, RoundDown)
    colors[indices]
end

function topo2net(topo::Topology)
    n, m = length(topo.nodes), length(topo.arcs)
    nodepos = reshape(reinterpret(Float64, topo.nodes), (2,n))'
    posx, posy = nodepos[:,1], nodepos[:,2]
    arcs = reshape(reinterpret(Int, topo.arcs), (2,m))
    Net(posx, posy, arcs)
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

    with(leg=false, grid=false, aspect_ratio=1,
         xlim=boundbox(posx), ylim=boundbox(posy)) do
        # draw edges
        plot(posx[edges], posy[edges], color=edgecolor, linewidth=edgewidth)

        # draw nodes
        scatter!(posx, posy, marker=(nodesize, nodecolor),
                 seriesann=map(x -> text(x, nodesize-3), 1:n))
    end
end

"drawing of networks, with mapping of attributes to color and size"
function draw(net::Net;
              edgecolor=:gray, edgecmap="Blues", edgewidth=4, edgebds=[3, 6],
              nodecolor=:orange, nodecmap="Reds", nodesize=12)
    if typeof(edgecolor) == Array{Float64, 2}
        edgecolor = map2color(edgecolor, cmap=edgecmap)
    end
    if typeof(edgewidth) == Array{Float64, 2}
        edgewidth = round(Int,normalize(edgewidth, lb=edgebds[1], ub=edgebds[2]))
    end
    if typeof(nodecolor) == Array{Float64, 2}
        nodecolor = map2color(nodecolor, cmap=nodecmap)
    end

    draw(net.posx, net.posy, net.edges,
         edgecolor=edgecolor, edgewidth=edgewidth,
         nodecolor=nodecolor, nodesize=nodesize)
end

draw(topo::Topology) = draw(topo2net(topo))

function draw(topo::Topology, cand::CandSol)
    ndiam = size(cand.zsol, 2)
    pipediam = cand.zsol * collect(linspace(1,ndiam,ndiam))

    actarcs = Arc[]
    actdiam = Float64[]
    for a = 1:length(topo.arcs)
        if pipediam[a] > 0
            push!(actarcs, topo.arcs[a])
            push!(actdiam, pipediam[a])
        end
    end

    net = topo2net(Topology(topo.nodes, actarcs))

    draw(net; edgewidth=actdiam', edgecolor=actdiam')
end

function draw(inst::Instance)
    net = topo2net(Topology(inst.nodes, []))

    px, py = net.posx, net.posy

    # area of nodes is relative to demand
    normdemand = inst.demand / maximum(abs(inst.demand))
    markersize = 20*[sqrt(abs(d)) for d in normdemand]
    markercolor = [d > 0 ? RGB(1,.5,.5) : RGB(.5,.5,1) for d in inst.demand]
    markerlabel = map(x -> text(x, 9), 1:length(inst.nodes))

    with(leg=false, grid=false, aspect_ratio=1,
         xlim=boundbox(px), ylim=boundbox(py)) do
        scatter(px, py, m=(markersize, markercolor), seriesann=markerlabel)
    end

    # TODO: draw pressure bounds
    # TODO: draw diameters
end

function draw_subsol(nodes::Vector{Node}, viol::Vector{Float64})
    n = length(nodes)
    nodepos = reshape(reinterpret(Float64, nodes), (2,n))'
    px, py = nodepos[:,1], nodepos[:,2]

    normviol = viol / maximum(abs(viol))
    markersize = 20*[sqrt(abs(v)) for v in normviol]
    markercolor = [v > 0 ? RGB(1,.5,.5) : RGB(.5,.5,1) for v in viol]
    markerlabel = map(x -> text(x, 9), 1:length(nodes))

    with(leg=false, grid=false, aspect_ratio=1,
         xlim=boundbox(px), ylim=boundbox(py)) do
        scatter(px, py, m=(markersize, markercolor), seriesann=markerlabel)
    end

    # TODO: maybe draw pressure loss on arcs?
end

function draw_dual(topo::Topology, cand::CandSol, λ::Vector{Float64},
                   μ::Vector{Float64})
    n,m = length(topo.nodes), length(topo.arcs)
    nodepos = reshape(reinterpret(Float64, topo.nodes), (2,n))'
    px, py = nodepos[:,1], nodepos[:,2]

    normλ = λ / maximum(abs(λ))
    markersize = 20*[sqrt(abs(v)) for v in normλ]
    markercolor = [v > 0 ? RGB(1,.5,.5) : RGB(.5,.5,1) for v in λ]
    markerlabel = map(x -> text(x, 9), 1:length(topo.nodes))

    candarcs = filter(a -> any(cand.zsol[a,:]), 1:m)
    μdense = fill(0, m)
    μdense[candarcs] = round(Int, 5 * μ / maximum(μ))
    arcs = reshape(reinterpret(Int, topo.arcs), (2,m))

    with(leg=false, grid=false, aspect_ratio=1,
         xlim=boundbox(px), ylim=boundbox(py)) do
        plot(px[arcs], py[arcs], linewidth=μdense', color="gray")
        scatter!(px, py, m=(markersize, markercolor), seriesann=markerlabel)
    end
end

function draw_paths(topo::Topology, paths, pathflows)
    nnodes, narcs = length(topo.nodes), length(topo.arcs)
    nodepos = reshape(reinterpret(Float64, topo.nodes), (2,nnodes))'
    px, py = nodepos[:,1], nodepos[:,2]

    # nodes
    sources = [path[1].tail for path in paths]
    sinks = [path[end].head for path in paths]
    @assert sources ∩ sinks == []

    markercolor = fill(RGB(1,1,1), nnodes)
    markercolor[sources] = RGB(.5,.5,1)
    markercolor[sinks] = RGB(1,.5,.5)

    markersize = fill(0.0, nnodes)
    markersize[sources] = 20
    markersize[sinks] = 20

    markerlabel = map(x -> text(x, 9), 1:nnodes)

    path2array(path) = reshape(reinterpret(Int, path), (2,length(path)))

    with(leg=false, grid=false, aspect_ratio=1, xlim=boundbox(px), ylim=boundbox(py)) do
        plot() # dummy initial plot
        for path in paths
            arcidx = path2array(path)
            plot!(px[arcidx], py[arcidx], color="gray", linewidth=4)
        end
        scatter!(px, py, m=(markersize, markercolor), seriesann=markerlabel)
    end
end
