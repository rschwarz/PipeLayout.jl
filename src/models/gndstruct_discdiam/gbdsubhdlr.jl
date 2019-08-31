"""
Constraint handler solving the GBD subproblem (fixed z).
"""
mutable struct GBDSubHdlr <: SCIP.AbstractConstraintHandler
    scip::SCIP.Optimizer # of master problem
    subsolver::JuMP.OptimizerFactory # for subproblems
    finaltime::Float64 # seconds
    addcritpath::Bool
    addnogoods::Bool
    counter::Int
    primal_bound::Float64
    best_solution::Union{CandSol, Nothing}

    function GBDSubHdlr(scip, subsolver, finaltime, addcritpaths, addnogoods)
        return new(scip, subsolver, finaltime, addcritpaths, addnogoods,
                   0, Inf, nothing)
    end
end

"""
Constraint data, referencing problem data and master variables.
"""
mutable struct GBDSubCons <: SCIP.AbstractConstraint{GBDSubHdlr}
    inst::Instance
    topo::Topology
    master::Master
end

function solve_gbd_sub(ch::GBDSubHdlr, cons::GBDSubCons,
                       sol::Ptr{SCIP.SCIP_SOL}=C_NULL; enforce=true)
    ch.counter += 1

    # fetch solution values
    ysol = SCIP.sol_values(ch.scip, cons.master.y, sol)
    zsol = SCIP.sol_values(ch.scip, cons.master.z, sol)
    qsol = SCIP.sol_values(ch.scip, cons.master.q, sol)
    ϕsol = SCIP.sol_values(ch.scip, cons.master.ϕ, sol)

    cand = CandSol(zsol .>= 0.5, qsol, ϕsol)

    # solve subproblem (from scratch, no warmstart)
    submodel, π, Δl, Δu, ploss, plb, pub =
        make_sub(cons.inst, cons.topo, cand, ch.subsolver)
    # TODO: solver.writemodels && writeLP(submodel, "sub_$(counter).lp", genericnames=false)
    settimelimit!(submodel, ch.subsolver, ch.finaltime - time())
    JuMP.optimize!(submodel)
    substatus = JuMP.termination_status(submodel)
    @assert substatus == MOI.OPTIMAL "Slack model is always feasible"
    totalslack = JuMP.objective_value(submodel)

    if totalslack ≈ 0.0
        # maybe only the relaxation is feasible, we have to check also the
        # "exact" subproblem with equations constraints.
        submodel2, _ = make_sub(cons.inst, cons.topo, cand, ch.subsolver,
                                relaxed=false)
        # TODO: solver.writemodels && writeLP(submodel2, "sub_exact_iter$(iter).lp", genericnames=false)
        settimelimit!(submodel2, ch.subsolver, ch.finaltime - time())
        JuMP.optimize!(submodel2)
        substatus2 = JuMP.termination_status(submodel2)
        @assert substatus2 == MOI.OPTIMAL "Slack model is always feasible"
        totalslack2 = JuMP.objective_value(submodel2)

        if totalslack2 ≈ 0.0
            # check whether this candidate is actually an improvement
            L = pipelengths(cons.topo)
            A = cons.topo.arcs
            D = cons.inst.diameters
            c = [diam.cost for diam in D]
            obj = [c[i] * L[a] for a=1:length(A), i=1:length(D)]
            candprimal = sum(obj .* zsol)
            if candprimal < ch.primal_bound
                # TODO: solver.debug && println("  found improving solution :-)")
                ch.best_solution = cand
                ch.primal_bound = candprimal
            else
                # TODO: solver.debug && println("  cand's obj val $(candprimal) no improvement over $(primal)!")
            end

            return SCIP.SCIP_FEASIBLE
        else
            # TODO: solver.debug && println("  subproblem/relaxation gap!")
            if !enforce
                return SCIP.SCIP_INFEASIBLE
            end

            # cut off candidate with no-good on z
            active = zsol .> 0.5
            coefs = 2.0 * active .- 1.0
            @cb_constraint(ch.scip, coefs ⋅ cons.master.z ≤ sum(active) - 1)

            return SCIP.SCIP_CONSADDED
        end
    else # totalslack > 0
        if !enforce
            return SCIP.SCIP_INFEASIBLE
        end

        # generate cuts and add to master
        if ch.addnogoods
            # cut off candidate with no-good on z
            active = zsol .> 0.5
            coefs = 2.0 * active .- 1.0
            @cb_constraint(ch.scip, coefs ⋅ cons.master.z ≤ sum(active) - 1)
        end

        if ch.addcritpath
            dualsol = SubDualSol(
                JuMP.dual.(ploss), JuMP.dual.(plb), JuMP.dual.(pub))
            critpathcuts(cons.inst, cons.topo, cons.master, cand,
                         subdualsol, ch.subsolver, cb=ch.scip)
        end

        # TODO: solver.debug && println("  added $(ncuts) cuts.")
        return SCIP.SCIP_CONSADDED
    end
end

function SCIP.check(ch::GBDSubHdlr,
                    constraints::Array{Ptr{SCIP.SCIP_CONS}},
                    sol::Ptr{SCIP.SCIP_SOL},
                    checkintegrality::Bool,
                    checklprows::Bool,
                    printreason::Bool,
                    completely::Bool)
    @assert length(constraints) == 1
    cons::GBDSubCons = SCIP.user_constraint(constraints[1])
    return solve_gbd_sub(ch, cons, sol, enforce=false)
end

function SCIP.enforce_lp_sol(ch::GBDSubHdlr,
                             constraints::Array{Ptr{SCIP.SCIP_CONS}},
                             nusefulconss::Cint,
                             solinfeasible::Bool)
    @assert length(constraints) == 1
    cons::GBDSubCons = SCIP.user_constraint(constraints[1])
    return solve_gbd_sub(ch, cons, enforce=true)
end

function SCIP.enforce_pseudo_sol(ch::GBDSubHdlr,
                                 constraints::Array{Ptr{SCIP.SCIP_CONS}},
                                 nusefulconss::Cint,
                                 solinfeasible::Bool,
                                 objinfeasible::Bool)
    @assert length(constraints) == 1
    cons::GBDSubCons = SCIP.user_constraint(constraints[1])
    return solve_gbd_sub(ch, cons, enforce=true)
end

function SCIP.lock(ch::GBDSubHdlr,
                   constraint::Ptr{SCIP.SCIP_CONS},
                   locktype::SCIP.SCIP_LOCKTYPE,
                   nlockspos::Cint,
                   nlocksneg::Cint)
    cons::GBDSubCons = SCIP.user_constraint(constraint)
    master::Master = cons.master
    for v in [master.y; master.z[:]; master.q; master.ϕ]
        # TODO: understand why lock is called during SCIPfree, after the
        # constraint should have been deleted already. Does it mean we should
        # implement CONSTRANS?
        var_::Ptr{SCIP.SCIP_VAR} = SCIP.var(ch.scip, JuMP.index(v))
        var_ != C_NULL || continue  # avoid segfault!

        SCIP.@SC SCIP.SCIPaddVarLocksType(
            ch.scip, var_, locktype, nlockspos + nlocksneg, nlockspos + nlocksneg)
    end
end
