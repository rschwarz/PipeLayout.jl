using JuMP

export ɛ, select_subset, settimelimit!, stilltime

# tolerance for numerical comparison
const ɛ = 1e-6

nonzero(x::Real) = abs(x) > ɛ

"select random subset of a base array of desired size (non-repeating)"
function select_subset{T}(base::Array{T}, size::Int)
    @assert size <= length(unique(base))
    candidates = T[]
    while length(candidates) < size
        candidates = unique(append!(candidates, rand(base, size)))
    end
    candidates[1:size]
end

const timebuffer = 10.0 # seconds

"update timelimit (in seconds) for internal solver"
function settimelimit!(model::JuMP.Model, limit)
    limit = max(limit, timebuffer) # at least buffer
    if limit < Inf
        setparameters!(internalmodel(model), TimeLimit=limit)
    end
end

"is there still enough until final limit?"
function stilltime(finaltime; buffer=timebuffer)
    time() + buffer < finaltime
end
