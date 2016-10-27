# utilities for network flow
using LightGraphs: is_connected, incidence_matrix, fadj, badj

export uniq_flow, flow_path_decomp

"Compute unique flow on arcs for given topology and demand vector."
function uniq_flow(topo::Topology, demand::Vector{Float64})
    length(topo.nodes) == length(demand) ||
        throw(ArgumentError("Topology and demand incompatible!"))
    is_tree(topo) || throw(ArgumentError("Topology is not a tree!"))
    sum(demand) ≈ 0.0 || throw(ArgumentError("Demand is not balanced!"))

    dg = digraph_from_topology(topo)
    N = incidence_matrix(dg, Float64)
    # TODO: figure out issues with CHOLMOD?
    full(N) \ demand
end

function uniq_flow(inst::Instance, topo::Topology)
    # allow for additional steiner nodes (demand = 0)
    termidx = [findfirst(topo.nodes, t) for t in inst.nodes]
    demand = fill(0.0, length(topo.nodes))
    demand[termidx] = inst.demand
    uniq_flow(topo, demand)
end

"Decompose arc flow into paths from sources to sinks."
function flow_path_decomp(topo::Topology, arcflow::Vector{Float64})
    dg = digraph_from_topology(topo)
    nnodes, narcs = length(topo.nodes), length(topo.arcs)
    @assert length(arcflow) == narcs
    @assert all(arcflow .≥ 0.0)

    arcidx = arcindex(topo)
    inarcs = [[arcidx[Arc(u,v)] for u in badj(dg,v)] for v in 1:nnodes]
    outarcs = [[arcidx[Arc(v,w)] for w in fadj(dg,v)] for v in 1:nnodes]

    # initialize residual flow
    resflow = copy(arcflow)
    resdem = fill(0.0, nnodes)
    for v = 1:nnodes
        resdem[v] += sum(resflow[inarcs[v]]) - sum(resflow[outarcs[v]])
    end
    @assert sum(resdem) ≈ 0.0

    paths = Vector{Arc}[]
    pathflows = Float64[]

    # iteratively move flow on paths
    while any(nonzero, resdem)
        path = []

        # start from largest remaining source
        current = indmin(resdem)
        pflow = -resdem[current]
        @assert pflow > 0.0

        while resdem[current] <= 0.0
            # pick largest outgoing arc
            maxarc = indmax([resflow[a] for a in outarcs[current]])
            curarc = outarcs[current][maxarc]
            pflow = min(pflow, resflow[curarc])
            @assert pflow > 0.0

            push!(path, topo.arcs[curarc])
            current = path[end].head
        end

        push!(paths, path)
        push!(pathflows, pflow)

        resdem[path[1].tail] += pflow
        for arc in path
            resflow[arcidx[arc]] -= pflow
        end
        resdem[path[end].head] -= pflow
    end

    @assert !any(nonzero, resdem)
    @assert !any(nonzero, resflow)

    paths, pathflows
end
