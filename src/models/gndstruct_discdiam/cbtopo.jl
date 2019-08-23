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
    ndiams = length(inst.diameters)
    narcs = length(topo.arcs)

    # reuse master problem from IterTopo algorithm
    model, y, q = make_semimaster(inst, topo, solver.mastersolver)
    solver.writemodels && writeLP(model, "master.lp", genericnames=false)

    # get access to SCIP.Optimizer instance (not factory)
    scip::SCIP.Optimizer = JuMP.backend(model)

    # enforce tree topology with constraint handler
    treetopohdlr = TreeTopoHdlr(scip)
    SCIP.include_conshdlr(scip, treetopohdlr; needs_constraints=true)
    SCIP.add_constraint(scip, treetopohdlr, TreeTopoCons(topo, y))

    # TODO: also add conshdlr for subproblem and other cuts
    # # callback counter and solution store
    # counter = 0
    # primal, bestsol = Inf, nothing

    # solve the master (including callback)
    JuMP.optimize!(model)
    status = JuMP.termination_status(model)
    if !(status in (MOI.OPTIMAL, MOI.INFEASIBLE, MOI.TIME_LIMIT, MOI.NODE_LIMIT))
        error("Unexpected status: $(status)")
    end

    # actually, we solved it, but the master does not know
    if status == MOI.INFEASIBLE && bestsol ≠ nothing
        status = MOI.OPTIMAL
    end
    solver.debug && println("  solved, status: $(status).")

    # TODO: extract dual bound, nnodes?
    Result(status, bestsol, primal, Inf, -1)
end
