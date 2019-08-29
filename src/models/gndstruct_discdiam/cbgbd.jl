"""
Callback-based algorithm for GBD in a ground structure, with discrete diameter
choice in master.

Solver object to store parameters.
"""
struct CallbackGBD <: GroundStructureSolver
    mastersolver
    subsolver
    timelimit::Float64 # seconds
    debug::Bool
    writemodels::Bool
    addcritpath::Bool
    addnogoods::Bool

    function CallbackGBD(mastersolver, subsolver;
                         timelimit=Inf, debug=false, writemodels=false,
                         addcritpath=true, addnogoods=false)
        new(mastersolver, subsolver, timelimit, debug, writemodels,
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
    scip::SCIP.Optimizer = JuMP.backend(model)

    # enforce tree topology with constraint handler (higher prio)
    treetopohdlr = TreeTopoHdlr(scip)
    SCIP.include_conshdlr(scip, treetopohdlr; needs_constraints=true,
                          enforce_priority=-15,
                          check_priority=-7_000_000)
    SCIP.add_constraint(scip, treetopohdlr, TreeTopoCons(topo, y))

    # enforce subproblem with constraint handler (lower prio)
    semisubhdlr = SemiSubHdlr(scip, solver.subsolver, finaltime)
    SCIP.include_conshdlr(
        scip, semisubhdlr;
        needs_constraints=true, enforce_priority=-16, check_priority=-7_100_000)
    SCIP.add_constraint(scip, semisubhdlr, SemiSubCons(
        inst, topo, master.y, master.z, master.q, master.ϕ))

    # solve master problem (has callback)
    JuMP.optimize!(master.model)
    status = JuMP.termination_status(master.model)

    # TODO: extract bounds, nnodes?
    dual = JuMP.objective_value(master.model)  # correct bound?
    nnodes = -1
    Result(status, bestsol, primal, dual, nnodes)
end
