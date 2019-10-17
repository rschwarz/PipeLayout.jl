"""
Callback-based algorithm for GBD in a ground structure, with discrete diameter
choice in master.

Solver object to store parameters.
"""
struct CallbackGBD <: GroundStructureSolver
    mastersolver
    subsolver
    timelimit::Float64 # seconds
    writemodels::Bool
    addcritpath::Bool
    addnogoods::Bool

    function CallbackGBD(mastersolver, subsolver;
                         timelimit=Inf, writemodels=false,
                         addcritpath=true, addnogoods=false)
        new(mastersolver, subsolver, timelimit, writemodels,
            addcritpath, addnogoods)
    end
end

"GBD on ground structure with callbacks"
function PipeLayout.optimize(inst::Instance, topo::Topology,
                             solver::CallbackGBD)
    # initialization
    finaltime = time() + solver.timelimit

    # no solution yet
    bestsol = nothing
    primal = Inf

    # use master problem from IterGBD algorithm
    master = Master(make_master(inst, topo, solver.mastersolver)...)
    solver.writemodels && writeLP(master.model, "master.lp", genericnames=false)

    # get access to SCIP.Optimizer instance (not factory)
    scip::SCIP.Optimizer = JuMP.backend(master.model)

    # enforce tree topology with constraint handler (higher prio)
    treetopohdlr = TreeTopoHdlr(scip)
    SCIP.include_conshdlr(
        scip, treetopohdlr;
        needs_constraints=true, enforce_priority=-15, check_priority=-7_000_000)
    SCIP.add_constraint(scip, treetopohdlr, TreeTopoCons(topo, master.y))

    # enforce subproblem with constraint handler (lower prio)
    gbdsubhdlr = GBDSubHdlr(scip, solver.subsolver, finaltime,
                            solver.addnogoods, solver.addcritpath)
    SCIP.include_conshdlr(
        scip, gbdsubhdlr;
        needs_constraints=true, enforce_priority=-16, check_priority=-7_100_000)
    SCIP.add_constraint(scip, gbdsubhdlr, GBDSubCons(inst, topo, master))

    # solve master problem (has callback)
    JuMP.optimize!(master.model)
    status = JuMP.termination_status(master.model)

    # TODO: correct bound?
    dual_bound = JuMP.objective_bound(master.model)
    # TODO: extract nnodes?
    nnodes = -1
    Result(status, gbdsubhdlr.best_solution, gbdsubhdlr.primal_bound,
           dual_bound, nnodes)
end
