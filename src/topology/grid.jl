export RegGrid, gridtopology

"Specification for M by N grid with node degree D"
immutable RegGrid{M,N,D}
    width::Float64
end

"Generate nodes of a regular grid"
function gridnodes{M,N,D}(grid::RegGrid{M,N,D})
    w = grid.width
    [Node(w*(j - 1), w*(M - i)) for j=1:N for i=1:M]
end

"2d array filled with indices to corresponding 1d array"
nodeidx(m,n) = reshape(1:m*n, m,n)

"generate coordinates for grid neighbors"
function candneighbors{M,N}(i::Int, j::Int, grid::RegGrid{M,N,4})
    [          (i-1, j),
     (i, j-1),           (i, j+1),
               (i+1, j)          ]
end
function candneighbors{M,N}(i::Int, j::Int, grid::RegGrid{M,N,8})
    [(i-1, j-1), (i-1, j), (i-1, j+1),
     (i, j-1),             (i, j+1),
     (i+1, j-1), (i+1, j), (i+1, j+1)]
end

"valid grid neighbors"
function gridneighbors{M,N,D}(i::Int, j::Int, grid::RegGrid{M,N,D})
    filter(t -> t[1] in 1:M && t[2] in 1:N,
           candneighbors(i, j, grid))
end

"Generate (anti-parallel) arcs of a regular grid"
function gridarcs{M,N,D}(grid::RegGrid{M,N,D})
    arcs = Arc[]
    idx = nodeidx(M,N)
    for j=1:N
        for i=1:M
            for n in gridneighbors(i, j, grid)
                push!(arcs, Arc(idx[i,j], idx[n...]))
            end
        end
    end
    arcs
end

function gridtopology{M,N,D}(grid::RegGrid{M,N,D})
    Topology(gridnodes(grid), gridarcs(grid))
end
