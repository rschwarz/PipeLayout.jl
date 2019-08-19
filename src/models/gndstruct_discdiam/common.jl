"Data type for candidate solutions from master problem."
struct CandSol
    zsol::Array{Bool,2}
    qsol::Vector{Float64}
    ϕsol::Vector{Float64}
end

"Result data from GBD algorithm"
struct Result
    status::Symbol
    solution
    primalbound::Float64
    dualbound::Float64
    niter::Int
end

"""
Extract sub-topology from active arcs and incident nodes.

The option `keep_nodes` allows to also include the isolated nodes.
"""
function topology_from_candsol(topo::Topology, ysol::Vector{Float64},
                               keep_nodes=false)
    actarcs = topo.arcs[ysol .> 0.5]
    if keep_nodes
        return Topology(topo.nodes, actarcs)
    end
    tails = [arc.tail for arc in actarcs]
    heads = [arc.head for arc in actarcs]
    nodeidx = tails ∪ heads
    nodemap = fill(0, length(topo.nodes))
    for (target, idx) in enumerate(nodeidx)
        nodemap[idx] = target
    end
    Topology(topo.nodes[nodeidx],
             [Arc(nodemap[arc.tail], nodemap[arc.head]) for arc in actarcs])
end

"Generic no-good cut"
function nogood(model::JuMP.Model, vars::AbstractArray{JuMP.VariableRef},
                sol::AbstractArray{T}) where T<:Real
    @assert size(vars) == size(sol)
    nvars = length(vars)
    active = (sol .> 0.5)
    nactive = sum(active)
    coef = 2.0*active .- 1.0
    @constraint(model, sum(coef[i]*vars[i] for i=1:nvars) <= nactive - 1)
    return 1
end

"Cut off all y values for given, undirected topology"
function avoid_topo_cut(model::JuMP.Model, y::AbstractArray{JuMP.VariableRef},
                        topo::Topology, edges::Vector{Arc})
    arcidx = arcindex(topo)
    antidx = antiparallelindex(topo)

    fwd = [arcidx[e] for e in edges]
    bwd = [antidx[a] for a in fwd]
    arcs = vcat(fwd, bwd)
    @constraint(model, sum(y[a] for a in arcs) ≤ length(edges) - 1)
end
