# custom types

immutable Node
    x::Float64 # [km]
    y::Float64 # [km]
end

immutable Bounds
    lb::Float64 # [bar]
    ub::Float64 # [bar]
end

immutable Diameter
    value::Float64 # [m]
    cost::Float64  # [EUR/m]
end

immutable Instance
    nodes::Vector{Node}
    demand::Vector{Float64} # [kg/s]
    pressure::Vector{Bounds}
    diameters::Vector{Diameter}
end

# JSON I/O
#  - serialization works out of the box,
#  - deserialization included:
include("deserialization.jl")

# Creating random instances
include("random.jl")
