"""
Constraint handler solving the semi-subproblem (fixed y, binary z).
"""
mutable struct SemiSubHdlr <: SCIP.AbstractConstraintHandler
    scip::SCIP.Optimizer
    solver::JuMP.OptimizerFactory
    finaltime::Float64 # seconds
    counter::Int
    primal_bound::Float64
    best_solution::Union{CandSol, Nothing}

    function SemiSubHdlr(scip, solver, finaltime)
        return new(scip, solver, finaltime, 0, Inf, nothing)
    end
end

"""
Constraint data, referencing the topology and binary edges variables.
"""
mutable struct SemiSubCons <: SCIP.AbstractConstraint{SemiSubHdlr}
    inst::Instance
    topo::Topology
    y::Array{JuMP.VariableRef}
    q::Array{JuMP.VariableRef}
end

function solve_semi_sub(ch::SemiSubHdlr, cons::SemiSubCons,
                        sol::Ptr{SCIP.SCIP_SOL}=C_NULL)
    ch.counter += 1

    # preparation
    inst = cons.inst
    topo = cons.topo
    ndiams = length(inst.diameters)
    narcs = length(topo.arcs)

    # fetch solution values
    ysol = SCIP.sol_values(ch.scip, cons.y, sol)
    qsol = SCIP.sol_values(ch.scip, cons.q, sol)
    # TODO: solver.debug && println("  cand. sol:$(find(qsol))")

    # complete solution candidate
    zcand = fill(false, narcs, ndiams)
    zcand[ysol .> 0.5, 1] .= true  # pretend to use smallest diameter
    cand = CandSol(zcand, qsol, qsol.^2)

    # solve subproblem (from scratch)
    submodel, candarcs, z = make_semisub(inst, topo, cand, ch.solver)
    # TODO: solver.writemodels && writeLP(submodel, "sub_$(counter).lp", genericnames=false)
    settimelimit!(submodel, ch.solver, ch.finaltime - time())
    JuMP.optimize!(submodel)
    substatus = JuMP.termination_status(submodel)
    if substatus == MOI.OPTIMAL
        # have found improving solution?
        newobj = JuMP.objective_value(submodel)
        if newobj < ch.primal_bound
            ch.primal_bound = newobj
            znew = fill(false, narcs, ndiams)
            znew[candarcs, :] = JuMP.value.(z) .> 0.5
            ch.best_solution = CandSol(znew, qsol, qsol.^2)
            # TODO: solver.debug && println("  new sol: $(primal)")
            return true # found improvement
        end
    elseif substatus != MOI.INFEASIBLE
        error("Unexpected status: $(:substatus)")
    end

    # TODO: what can we do about dual/primal stopping criterion?
    # TODO: add cut on objective var?

    return false # did not find improvement
end

function SCIP.check(ch::SemiSubHdlr,
                    constraints::Array{Ptr{SCIP.SCIP_CONS}},
                    sol::Ptr{SCIP.SCIP_SOL},
                    checkintegrality::Bool,
                    checklprows::Bool,
                    printreason::Bool,
                    completely::Bool)
    @assert length(constraints) == 1
    cons::SemiSubCons = SCIP.user_constraint(constraints[1])

    # Always solve semi-subproblem (to find improving solutions)
    # TODO: cache candidates to avoid re-solving?
    improvement = solve_semi_sub(ch, cons, sol)

    # Never accept the solution candidate.
    return SCIP.SCIP_INFEASIBLE
end


function enforce_semi_sub(ch::SemiSubHdlr, cons::SemiSubCons)
    # Always solve semi-subproblem (to find improving solutions)
    # TODO: cache candidates to avoid re-solving?
    improvement = solve_semi_sub(ch, cons)

    # Always add no-good to cut off this candidate.
    ysol = SCIP.sol_values(ch.scip, cons.y)
    active = ysol .> 0.5
    coefs = 2.0 * active .- 1.0
    ci = add_cons(ch.scip, @build_constraint(
        coefs ⋅ cons.y ≤ sum(active) - 1))

    return SCIP.SCIP_CONSADDED
end

function SCIP.enforce_lp_sol(ch::SemiSubHdlr,
                             constraints::Array{Ptr{SCIP.SCIP_CONS}},
                             nusefulconss::Cint,
                             solinfeasible::Bool)
    @assert length(constraints) == 1
    cons::SemiSubCons = SCIP.user_constraint(constraints[1])
    return enforce_semi_sub(ch, cons)
end

function SCIP.enforce_pseudo_sol(ch::SemiSubHdlr,
                                 constraints::Array{Ptr{SCIP.SCIP_CONS}},
                                 nusefulconss::Cint,
                                 solinfeasible::Bool,
                                 objinfeasible::Bool)
    @assert length(constraints) == 1
    cons::SemiSubCons = SCIP.user_constraint(constraints[1])
    return enforce_semi_sub(ch, cons)
end

function SCIP.lock(ch::SemiSubHdlr,
                   constraint::Ptr{SCIP.SCIP_CONS},
                   locktype::SCIP.SCIP_LOCKTYPE,
                   nlockspos::Cint,
                   nlocksneg::Cint)
    cons::SemiSubCons = SCIP.user_constraint(constraint)
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
