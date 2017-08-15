"""
Callback-based algorithm for decomposition on topology for Ground Structure with
Discrete Diameters.

Solver object to store parameter values.
"""
struct CallbackTopo <: GroundStructureSolver
    mastersolver
    subsolver
    timelimit::Float64 # seconds
    debug::Bool
    writemodels::Bool

    function CallbackTopo(mastersolver, subsolver;
                          timelimit=Inf, debug=False, writemodels=False)
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
    writemodels && writeLP(model, "master.lp")

    # callback counter and solution store
    counter = 0
    primal, bestsol = Inf, nothing

    # add callback that solves subproblem
    function cbtopo(cb)
        ysol, qsol = getvalue(y), getvalue(q)
        debug && println("  cand. sol:$(find(qsol))")
        counter += 1

        # TODO: reduce code duplication with itertopo.jl
        # check whether candidate has tree topology
        candtopo = topology_from_candsol(topo, ysol)
        if !is_tree(candtopo)
            fullcandtopo = topology_from_candsol(topo, ysol, true)
            cycle = find_cycle(fullcandtopo)
            if length(cycle) == 0
                # TODO: deal with disconnected candidate?
                debug && println("  cut off disconnected topology with nogood.")
                nogood(mastermodel, y, ysol, cb=cb)
            else
                debug && println("  cut off non-tree topology, cycle: $(cycle)")
                avoid_topo_cut(mastermodel, y, topo, cycle, cb=cb)
            end
            return  # exit callback
        end

        # solve subproblem (from scratch, no warmstart)
        zcand = fill(false, narcs, ndiams)
        zcand[ysol .> 0.5, 1] = true
        cand = CandSol(zcand, qsol, qsol.^2)

        submodel, candarcs, z = make_semisub(inst, topo, cand, solver.subsolver)
        writemodels && writeLP(submodel, "sub_$(counter).lp")
        settimelimit!(submodel, solver.subsolver, finaltime - time())
        substatus = solve(submodel, suppress_warnings=true)
        if substatus == :Optimal
            # have found improving solution?
            newobj = getobjectivevalue(submodel)
            if newobj < primal
                primal = newobj
                znew = fill(false, narcs, ndiams)
                znew[candarcs, :] = (getvalue(z) .> 0.5)
                bestsol = CandSol(znew, qsol, qsol.^2)
                debug && println("  found improving solution: $(primal)")
            end
        elseif substatus != :Infeasible
            error("Unexpected status: $(:substatus)")
        end

        # TODO: what can we do about dual/primal stopping criterion?
        # TODO: add cut on objective var?

        # generate nogood cut and add to master
        nogood(model, y, ysol, cb)
    end
    addlazycallback(model, cbtopo)

    # solve the master (including callback)
    status = solve(mastermodel, suppress_warnings=true)
    if status ≠ :Optimal && status ≠ :Infeasible && status ≠ :UserLimit
        error("Unexpected status: $(status)")
    end

    # actually, we solved it, but the master does not know
    if status == :Infeasible && bestsol ≠ nothing
        status = :Optimal
    end
    debug && println("  problem solved with status $(status).")

    # TODO: extract dual bound, nnodes?
    Result(status, bestsol, primal, Inf, nothing)
end
