export ɛ, select_subset

# tolerance for numerical comparison
const ɛ = 1e-6

nonzero(x::Real) = abs(x) > ɛ

" select random subset of a base array of desired size (non-repeating)"
function select_subset{T}(base::Array{T}, size::Int)
    @assert size <= length(unique(base))
    candidates = T[]
    while length(candidates) < size
        candidates = unique(append!(candidates, rand(base, size)))
    end
    candidates[1:size]
end
