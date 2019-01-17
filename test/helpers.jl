#
# test data
#
"create a small instance (and topology) with a single pipe"
function single_pipe(;length=100.0, flow=200.0)
    nodes = [Node(0,0), Node(length,0)]
    demand = [-flow, flow]
    press = fill(Bounds(60.0, 80.0), size(nodes))
    diams = [Diameter(t...) for t in [(0.8, 1.0),(1.0, 1.2),(1.2, 1.5)]]
    inst = Instance(nodes, demand, press, diams)
    topo = Topology(nodes, [Arc(1,2)])
    inst, topo
end

"Create a grid of square cells with m by n nodes and given edge width."
function squaregrid(m::Int, n::Int, width::T; antiparallel=false) where T<:Real
    nodes = Node[]
    sizehint!(nodes, m*n)
    for j in 1:n
        x = width * (j - 1)
        for i in 1:m
            push!(nodes, Node(x, width * (m - i)))
        end
    end
    nodeidx = reshape(1:m*n, (m,n))

    arcs = Arc[]
    narcs = (m-1)*n + m*(n-1)
    sizehint!(arcs, antiparallel ? 2*narcs : narcs)
    # columns top-down
    for j in 1:n
        for i in 1:m-1
            push!(arcs, Arc(nodeidx[i,j], nodeidx[i+1,j]))
            antiparallel && push!(arcs, Arc(nodeidx[i+1,j], nodeidx[i,j]))
        end
    end
    # rows left-right
    for i in 1:m
        for j in 1:n-1
            push!(arcs, Arc(nodeidx[i,j], nodeidx[i,j+1]))
            antiparallel && push!(arcs, Arc(nodeidx[i,j+1], nodeidx[i,j]))
        end
    end

    Topology(nodes, arcs)
end
