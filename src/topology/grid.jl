export squaregrid

"Create a grid of square cells with m by n nodes and given edge width."
function squaregrid{T<:Real}(m::Int, n::Int, width::T; antiparallel=false)
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
