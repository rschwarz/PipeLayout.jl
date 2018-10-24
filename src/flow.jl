# utilities for network flow
using LightGraphs: inneighbors, outneighbors

export uniq_flow, flow_path_decomp

"Create dense demand vector for a topology"
function dense_demand(inst::Instance, topo::Topology)
    if length(inst.nodes) == length(topo.nodes)
        # just stop if there are no extra Steiner nodes
        return inst.demand
    end
    if sum(topo.nodes .== topo.nodes') > length(topo.nodes)
        # there must be duplicate nodes, leading to ambiguity!
        throw(ArgumentError("Duplicate nodes in topology!"))
    end

    # allow for additional steiner nodes (demand = 0)
    termidx = termindex(topo.nodes, inst.nodes)
    demand = fill(0.0, length(topo.nodes))
    demand[termidx] = inst.demand
    demand
end

"Compute unique flow on arcs for given topology and demand vector."
function uniq_flow(topo::Topology, demand::Vector{Float64})
    length(topo.nodes) == length(demand) ||
        throw(ArgumentError("Topology and demand incompatible!"))
    is_tree(topo) || throw(ArgumentError("Topology is not a tree!"))
    sum(demand) ≈ 0.0 || throw(ArgumentError("Demand is not balanced!"))

    N = incidence(topo)
    N \ demand
end

function uniq_flow(inst::Instance, topo::Topology)
    uniq_flow(topo, dense_demand(inst, topo))
end

"Reorient the arcs of a tree so that flow goes forward (positive)"
function reorient_fwdflow(topo::Topology, demand::Vector{Float64})
    arcflow = uniq_flow(topo, demand)
    arcs = Arc[]
    for (arc, flow) in zip(topo.arcs, arcflow)
        if flow ≥ 0
            push!(arcs, arc)
        else
            push!(arcs, Arc(arc.head, arc.tail))
        end
    end
    Topology(topo.nodes, arcs)
end

function reorient_fwdflow(inst::Instance, topo::Topology)
    reorient_fwdflow(topo, dense_demand(inst, topo))
end

"Decompose arc flow into paths from sources to sinks."
function flow_path_decomp(topo::Topology, arcflow::Vector{Float64})
    dg = digraph_from_topology(topo)
    nnodes, narcs = length(topo.nodes), length(topo.arcs)
    @assert length(arcflow) == narcs
    @assert all(arcflow .≥ 0.0)

    arcidx = arcindex(topo)
    inarcs = [[arcidx[Arc(u,v)] for u in inneighbors(dg,v)] for v in 1:nnodes]
    outarcs = [[arcidx[Arc(v,w)] for w in outneighbors(dg,v)] for v in 1:nnodes]

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
        current = argmin(resdem)
        pflow = -resdem[current]
        @assert pflow > 0.0

        while resdem[current] <= 0.0
            # pick largest outgoing arc
            maxarc = argmax([resflow[a] for a in outarcs[current]])
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
