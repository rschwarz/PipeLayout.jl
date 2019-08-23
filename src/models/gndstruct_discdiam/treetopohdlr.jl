"""
Constraint handler enforcing tree topologies.
"""
mutable struct TreeTopoHdlr <: SCIP.AbstractConstraintHandler
    scip::SCIP.Optimizer
end

"""
Constraint data, referencing the topology and binary edges variables.
"""
mutable struct TreeTopoCons <: SCIP.AbstractConstraint{TreeTopoHdlr}
    topo::Topology
    y::Array{JuMP.VariableRef}
end

function SCIP.check(ch::TreeTopoHdlr,
                    constraints::Array{Ptr{SCIP.SCIP_CONS}},
                    sol::Ptr{SCIP.SCIP_SOL},
                    checkintegrality::Bool,
                    checklprows::Bool,
                    printreason::Bool,
                    completely::Bool)
    @assert length(constraints) == 1
    cons::TreeTopoCons = SCIP.user_constraint(constraints[1])
    ysol = SCIP.sol_values(ch.scip, cons.y, sol)

    # no need to find out more details
    return (is_tree(topology_from_candsol(cons.topo, ysol))
            ? SCIP.SCIP_FEASIBLE
            : SCIP.SCIP_INFEASIBLE)
end

"""
Enforce tree topology.

Case 1: Topology is disconnected --> add no-good on all edges.
Case 2: Topology contains cycle --> add no-good on cycle edges.
"""
function enforce_tree_topo(ch::TreeTopoHdlr, cons::TreeTopoCons)
    ysol = SCIP.sol_values(ch.scip, cons.y)

    # Tree topology is fine, nothing to do.
    if is_tree(topology_from_candsol(cons.topo, ysol))
        return SCIP.SCIP_FEASIBLE
    end

    # Find cycle in topology.
    cycle = find_cycle(topology_from_candsol(cons.topo, ysol, true))

    # Case distinction
    if length(cycle) == 0
        # 1) There is no cycle. So, topology must be disconnected (forest).
        # TODO: Refactor no-good specific code for reuse.
        # TODO: Can we do better than no-good by forcing cut edges only?
        active = ysol .> 0.5
        coefs = 2.0 * active .- 1.0
        cut = @build_constraint(coefs ⋅ cons.y ≤ sum(active) - 1)
        ci = MOI.add_constraint(ch.scip, cut)
    else
        # 2) There is a cycle. Let's forbid it (with anti-parallel arcs).
        arcidx = arcindex(cons.topo)
        antidx = antiparallelindex(cons.topo)
        fwd = [arcidx[e] for e in cycle]
        bwd = [antiparallelindex[e] for e in fwd]
        arcs = vcat(fwd, bwd)
        cut = @build_constraint(sum(cons.y[a] for a in arcs) ≤ length(cycle) - 1)
        ci = MOI.add_constraint(ch.scip, cut)
    end

    return SCIP.SCIP_CONSADDED
end

function SCIP.enforce_lp_sol(ch::TreeTopoHdlr,
                             constraints::Array{Ptr{SCIP.SCIP_CONS}},
                             nusefulconss::Cint,
                             solinfeasible::Bool)
    @assert length(constraints) == 1
    cons::TreeTopoCons = SCIP.user_constraint(constraints[1])
    return enforce_tree_topo(ch, cons)
end

function SCIP.enforce_pseudo_sol(ch::TreeTopoHdlr,
                                 constraints::Array{Ptr{SCIP.SCIP_CONS}},
                                 nusefulconss::Cint,
                                 solinfeasible::Bool,
                                 objinfeasible::Bool)
    @assert length(constraints) == 1
    cons::TreeTopoCons = SCIP.user_constraint(constraints[1])
    return enforce_tree_topo(ch, cons)
end

function SCIP.lock(ch::TreeTopoHdlr,
                   constraint::Ptr{SCIP.SCIP_CONS},
                   locktype::SCIP.SCIP_LOCKTYPE,
                   nlockspos::Cint,
                   nlocksneg::Cint)
    cons::TreeTopoCons = SCIP.user_constraint(constraint)
    for y in cons.y
        # TODO: understand why lock is called during SCIPfree, after the
        # constraint should have been deleted already. Does it mean we should
        # implement CONSTRANS?
        var_::Ptr{SCIP.SCIP_VAR} = SCIP.var(ch.scip, JuMP.index(y))
        var_ != C_NULL || continue  # avoid segfault!

        SCIP.@SC SCIP.SCIPaddVarLocksType(
            ch.scip, var_, locktype, nlockspos + nlocksneg, nlockspos + nlocksneg)
    end
end
