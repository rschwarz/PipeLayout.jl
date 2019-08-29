"""
Constraint handler solving the GBD subproblem (fixed z).
"""
mutable struct GBDSubHdlr <: SCIP.AbstractConstraintHandler
    scip::SCIP.Optimizer # of master problem
    solver::JuMP.OptimizerFactory # for subproblems
    finaltime::Float64 # seconds
    counter::Int

    function GBDSubHdlr(scip, solver, finaltime)
        return new(scip, solver, finaltime, 0)
    end
end

"""
Constraint data, referencing problem data and master variables.
"""
mutable struct GBDSubCons <: AbstractConstraint{GBDSubHdlr}
    inst::Instance
    topo::Topology
    y::Array{JuMP.VariableRef} # arc selection
    z::Array{JuMP.VariableRef} # diameter selection
    q::Array{JuMP.VariableRef} # flow on arc
end

function solve_gbd_sub(ch::GBDSubHdlr, cons::GBDSubCons,
                       sol::Ptr{SCIP.SCIP_SOL}=C_NULL, enforce=true)
    ch.counter += 1

    # preparation

    # fetch solution values
    ysol = SCIP.sol_values(ch.scip, cons.y, sol)
    zsol = SCIP.sol_values(ch.scip, cons.y, sol)
    qsol = SCIP.sol_values(ch.scip, cons.q, sol)
    ϕsol = SCIP.sol_values(ch.scip, cons.ϕ, sol)

    cand = CandSol(zcand .>= 0.5, qsol, ϕsol)

    # solve subproblem (from scratch, no warmstart)
    submodel, π, Δl, Δu, ploss, plb, pub =
        make_sub(cons.inst, cons.topo, cand, ch.solver)
    solver.writemodels && writeLP(submodel, "sub_$(counter).lp", genericnames=false)
    settimelimit!(submodel, subsolver, ch.finaltime - time())
    JuMP.optimize!(submodel)
    substatus = JuMP.termination_status(submodel)
    @assert substatus == MOI.OPTIMAL "Slack model is always feasible"
    totalslack = JuMP.objective_value(submodel)
    if totalslack ≈ 0.0
        # maybe only the relaxation is feasible, we have to check also the
        # "exact" subproblem with equations constraints.
        submodel2, _ = make_sub(cons.inst, cons.topo, cand, ch.solver, relaxed=false)
        solver.writemodels && writeLP(submodel2, "sub_exact_iter$(iter).lp", genericnames=false)
        settimelimit!(submodel2, ch.solver, ch.finaltime - time())
        JuMP.optimize!(submodel2)
        substatus2 = JuMP.termination_status(submodel2)
        @assert substatus2 == MOI.OPTIMAL "Slack model is always feasible"
        totalslack2 = JuMP.objective_value(submodel2)

        if totalslack2 ≈ 0.0
            # check whether this candidate is actually an improvement
            L = pipelengths(topo)
            c = [diam.cost for diam in inst.diameters]
            obj = [c[i] * L[a] for a=1:length(topo.arcs), i=1:length(inst.diameters)]
            candprimal = sum(obj .* z)
            if candprimal < primal
                solver.debug && println("  found improving solution :-)")
                bestsol = cand
                primal = candprimal
            elseif solver.debug
                println("  cand's obj val $(candprimal) no improvement over $(primal)!")
            end
        else
            # cut off candidate with no-good on z
            solver.debug && println("  subproblem/relaxation gap!")
            nogood(master.model, master.z, cand.zsol)
        end

        return  # exit callback
    end

    dualsol = SubDualSol(JuMP.dual.(ploss), JuMP.dual.(plb), JuMP.dual.(pub))

    # generate cuts and add to master
    ncuts = cuts(inst, topo, master, cand, dualsol, solver.subsolver,
                 addnogoods=solver.addnogoods,
                 addcritpath=solver.addcritpath, cb=cb)
    solver.debug && println("  added $(ncuts) cuts.")


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

end

function enforce_gbd_sub(ch::GBDSubHdlr, cons::GBDSubCons)


    return SCIP.SCIP_CONSADDED
end

function SCIP.enforce_lp_sol(ch::GBDSubHdlr,
                             constraints::Array{Ptr{SCIP.SCIP_CONS}},
                             nusefulconss::Cint,
                             solinfeasible::Bool)
    @assert length(constraints) == 1
    cons::GBDSubCons = SCIP.user_constraint(constraints[1])
    return enforce_gbd_sub(ch, cons)
end

function SCIP.enforce_pseudo_sol(ch::GBDSubHdlr,
                                 constraints::Array{Ptr{SCIP.SCIP_CONS}},
                                 nusefulconss::Cint,
                                 solinfeasible::Bool,
                                 objinfeasible::Bool)
    @assert length(constraints) == 1
    cons::GBDSubCons = SCIP.user_constraint(constraints[1])
    return enforce_gbd_sub(ch, cons)
end

function SCIP.lock(ch::GBDSubHdlr,
                   constraint::Ptr{SCIP.SCIP_CONS},
                   locktype::SCIP.SCIP_LOCKTYPE,
                   nlockspos::Cint,
                   nlocksneg::Cint)
    cons::GBDSubCons = SCIP.user_constraint(constraint)
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
