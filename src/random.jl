"""Sample (uniformly) random node locations within [0, xmax]*[0, ymax].
Round coordinates to integer multiples of grid."""
function randnodes(num::Int64, xmax::Int64, ymax::Int64, grid::Int64=1)
    xcoords = rand(0:grid:xmax, num)
    ycoords = rand(0:grid:ymax, num)
    [Node(t...) for t in zip(xcoords, ycoords)]
end

"random probability distribution on n events"
function randdist(n::Int64)
    x = rand(n) + 0.1*ones(n)
    x / norm(x, 1)
end

"""random demand vector for n nodes, with indices of entries and
exits, and total demand"""
function randdemand(n::Int64, entries::Vector{Int64}, exits::Vector{Int64},
                    totalflow::Float64)
    demand = zeros(n)
    demand[entries] = - totalflow * randdist(length(entries))
    demand[exits] = totalflow * randdist(length(exits))
    demand
end

"Randomly partition points into entries and exits, with at least one of each"
function randboundary(num::Int64)
    @assert num >= 2
    indicator = rand(Bool, num - 2)
    rest = collect(3:num)
    entries = vcat([1], rest[indicator])
    exits = vcat([2], rest[!indicator])
    entries, exits
end


"""random instance on a grid"""
function randgridinstance(num::Int, xmax::Int64, ymax::Int64, flow::Float64,
                          diameters::Vector{Diameter}, grid::Int64=5,
                          lb::Float64=40.0, ub::Float64=90.0)
    entries, exits = randboundary(num)
    Instance(randnodes(num, xmax, ymax, grid),
             randdemand(num, entries, exits, flow),
             fill(Bounds(lb, ub), num),
             diameters)
end

"select random subset of a base array of desired size (non-repeating)"
function select_subset{T}(base::Array{T}, size::Int)
    @assert size <= length(unique(base))
    candidates = T[]
    while length(candidates) < size
        candidates = unique(append!(candidates, rand(base, size)))
    end
    candidates[1:size]
end

"rejection sampling for good area coverage"
function select_covering_subset(base::Array{Node}, size::Int)
    function area(nodes::Array{Node})
        xl, xu = extrema([n.x for n in nodes])
        yl, yu = extrema([n.y for n in nodes])
        (xu - xl)*(yu - yl)
    end

    whole = area(base)
    fraction = size/(size+1)
    sample = base[1:1]
    while area(sample) < fraction * whole
        sample = select_subset(base, size)
    end
    sample
end
