# custom types

export Node, Bounds, Diameter, Instance, Arc, Topology

"Node type with location on plane."
immutable Node
    x::Float64 # [km]
    y::Float64 # [km]
end

"Variable bounds as interval, used for node pressure."
immutable Bounds
    lb::Float64 # [bar]
    ub::Float64 # [bar]
end

"Diameter has value and cost factor (per length)."
immutable Diameter
    value::Float64 # [m]
    cost::Float64  # [EUR/m]
end

"""
An instance for the Pipe Layout Problem.

It is defined by boundary nodes, a fixed and balanced flow demand, pressure
bounds and available diameters.
"""
immutable Instance
    nodes::Vector{Node}
    demand::Vector{Float64} # [kg/s]
    pressure::Vector{Bounds}
    diameters::Vector{Diameter}
end

"An arc is specified by integer indices to the node array"
immutable Arc
    tail::Int
    head::Int
end

"""
A topology for a pipe layout.

May contain (Steiner) nodes (in addition to the boundary nodes defined in the
instance). These would be appended at the end. Their location need not be
considered fixed.

Contains arcs that are given as pairs of indices into the array of nodes.

Can be used both as a candidate topology (in form of a spanning tree) or as the
ground structure for the choice of topologies.
"""
immutable Topology
    nodes::Vector{Node}
    arcs::Vector{Arc}
end

# JSON I/O
#  - serialization works out of the box,
#  - deserialization included:
include("deserialization.jl")

# Creating random instances
include("random.jl")
