using JuMP
using MathOptInterface

const MOI = MathOptInterface

struct Solution
    nodes::Vector{Node}
    lsol::Array{Float64,2}
    πsol::Vector{Float64}
end

struct Result
    status::MOI.TerminationStatusCode
    sol::Solution
    value::Float64
end

include("nlp.jl")
include("soc.jl")
