# utilities for network flow
using LightGraphs: is_connected, incidence_matrix

"Compute unique flow on arcs for given topology and demand vector."
function uniq_flow(topo::Topology, demand::Vector{Float64})
    length(topo.nodes) == length(demand) ||
        throw(ArgumentError("Topology and demand incompatible!"))
    is_tree(topo) || throw(ArgumentError("Topology is not a tree!"))
    sum(demand) ≈ 0.0 || throw(ArgumentError("Demand is not balanced!"))

    dg = digraph_from_topology(topo)
    incidence_matrix(dg) \ demand
end

uniq_flow(inst::Instance, topo::Topology) = uniq_flow(topo, inst.demand)
