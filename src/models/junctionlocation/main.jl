using JuMP
using MathProgBase

immutable Solution
    nodes::Vector{Node}
    lsol::Array{Float64,2}
    πsol::Vector{Float64}
end

immutable Result
    status::Symbol
    sol::Solution
    value::Float64
end

include("nlp.jl")
include("soc.jl")
