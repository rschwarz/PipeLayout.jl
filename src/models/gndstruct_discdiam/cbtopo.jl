"""
Callback-based algorithm for decomposition on topology for Ground Structure with
Discrete Diameters.

Solver object to store parameter values.
"""
struct CallbackTopo <: GroundStructureSolver
    mastersolver::JuMP.OptimizerFactory
    subsolver::JuMP.OptimizerFactory
    timelimit::Float64 # seconds
    debug::Bool
    writemodels::Bool

    function CallbackTopo(mastersolver, subsolver;
                          timelimit=Inf, debug=false, writemodels=false)
        new(mastersolver, subsolver, timelimit, debug, writemodels)
    end
end

"Callbacks for decomposition based on topology"
function PipeLayout.optimize(inst::Instance, topo::Topology,
                             solver::CallbackTopo)
    # initialization
    finaltime = time() + solver.timelimit

    # reuse master problem from IterTopo algorithm
    model, y, q = make_semimaster(inst, topo, solver.mastersolver)
    solver.writemodels && writeLP(model, "master.lp", genericnames=false)

    # get access to SCIP.Optimizer instance (not factory)
    scip::SCIP.Optimizer = JuMP.backend(model)

    # enforce tree topology with constraint handler (higher prio)
    treetopohdlr = TreeTopoHdlr(scip)
    SCIP.include_conshdlr(scip, treetopohdlr; needs_constraints=true,
                          enforce_priority=-15,
                          check_priority=-7_000_000)
    SCIP.add_constraint(scip, treetopohdlr, TreeTopoCons(topo, y))

    # enforce feasible subproblems with constraint handler (lower prio)
    semisubhdlr = SemiSubHdlr(scip, solver.subsolver, finaltime)
    SCIP.include_conshdlr(scip, semisubhdlr; needs_constraints=true,
                          enforce_priority=-16,
                          check_priority=-7_100_000)
    SCIP.add_constraint(scip, semisubhdlr,
                        SemiSubCons(inst, topo, y, q))

    # solve the master (including callback)
    JuMP.optimize!(model)
    status = JuMP.termination_status(model)
    if !(status in (MOI.OPTIMAL, MOI.INFEASIBLE, MOI.TIME_LIMIT, MOI.NODE_LIMIT))
        error("Unexpected status: $(status)")
    end

    # actually, we solved it, but the master does not know
    if status == MOI.INFEASIBLE && semisubhdlr.best_solution !== nothing
        status = MOI.OPTIMAL
    end
    solver.debug && println("  solved, status: $(status).")

    # TODO: extract dual bound, nnodes?
    Result(status, semisubhdlr.best_solution, semisubhdlr.primal_bound, Inf, -1)
end
