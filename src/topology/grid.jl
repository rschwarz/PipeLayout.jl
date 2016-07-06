"Create a grid of square cells with m by n nodes and given edge width."
function squaregrid(m::Int, n::Int, width::Float64)
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
    sizehint!(arcs, (m-1)*n + m*(n-1))
    # columns top-down
    for j in 1:n
        for i in 1:m-1
            push!(arcs, [nodeidx[i,j], nodeidx[i+1,j]])
        end
    end
    # rows left-right
    for i in 1:m
        for j in 1:n-1
            push!(arcs, [nodeidx[i,j], nodeidx[i,j+1]])
        end
    end

    Topology(nodes, arcs)
end
