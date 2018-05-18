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
    counter = 0
    bestsol = nothing

    # use master problem from IterGBD algorithm
    master = Master(make_master(inst, topo, solver.mastersolver)...)
    solver.writemodels && writeLP(master.model, "master.lp", genericnames=false)

    # TODO: rm code duplication with itergbd
    # add callback for subproblem & cuts
    function cbgbd(cb)
        ysol = getvalue(master.y)
        zsol = getvalue(master.z)
        cand = CandSol(zsol .>= 0.5,
                       getvalue(master.q),
                       getvalue(master.ϕ))
        if solver.debug
            j,i,_ = findnz(cand.zsol')
            println("  cand. sol:$(collect(zip(i,j)))")
        end
        counter += 1

        # check whether candidate has tree topology
        candtopo = topology_from_candsol(topo, ysol)
        if !is_tree(candtopo)
            fullcandtopo = topology_from_candsol(topo, ysol, true)
            cycle = find_cycle(fullcandtopo)
            if length(cycle) == 0
                # TODO: Actually, this might be optimal, but it could also occur
                # when adding some irrelevant pipe is cheaper than increasing
                # the diameter. How to distinguish these cases?
                nogood(master.model, master.y, ysol, cb=cb)
                solver.debug && println("  cb: nogood cut: discon topology.")
            else
                avoid_topo_cut(master.model, master.y, topo, cycle, cb=cb)
                solver.debug && println("  cb: cut topology w/ cycle: $(cycle)")
            end
            return  # exit from callback
        end

        # check whether y and z are consistent
        zsum = sum(zsol, 2)
        yz_match = (zsum .> 0.5) .== (ysol .> 0.5)
        if !all(yz_match)
            if solver.debug
                i, _, _ = findnz(.!yz_match)
                println("  cb: mismatch between y and z values for arcs $(i)")
                println("    ysol = $(ysol[i])")
                println("    zsol = $(zsol[i,:])")
            end
            return
        end

        # solve subproblem (from scratch, no warmstart)
        submodel, π, Δl, Δu, ploss, plb, pub =
            make_sub(inst, topo, cand, solver.subsolver)
        solver.writemodels && writeLP(submodel, "sub_$(counter).lp", genericnames=false)
        settimelimit!(submodel, solver.subsolver, finaltime - time())
        substatus = solve(submodel, suppress_warnings=true)
        @assert substatus == :Optimal "Slack model is always feasible"
        totalslack = getobjectivevalue(submodel)
        if totalslack ≈ 0.0
            # maybe only the relaxation is feasible, we have to check also the
            # "exact" subproblem with equations constraints.
            submodel2, _ = make_sub(inst, topo, cand, solver.subsolver, relaxed=false)
            solver.writemodels && writeLP(submodel2, "sub_exact_iter$(iter).lp", genericnames=false)
            settimelimit!(submodel2, solver.subsolver, finaltime - time())
            substatus2 = solve(submodel2, suppress_warnings=true)
            @assert substatus2 == :Optimal "Slack model is always feasible"
            totalslack2 = getobjectivevalue(submodel2)

            if totalslack2 ≈ 0.0
                solver.debug && println("  found feasible solution :-)")
                bestsol = cand
                # don't need to do anything, solved the problem
            else
                # cut off candidate with no-good on z
                solver.debug && println("  subproblem/relaxation gap!")
                nogood(master.model, master.z, cand.zsol)
            end

            return  # exit callback
        end

        dualsol = SubDualSol(getdual(ploss), getdual(plb), getdual(pub))

        # generate cuts and add to master
        ncuts = cuts(inst, topo, master, cand, dualsol, solver.subsolver,
                     addnogoods=solver.addnogoods,
                     addcritpath=solver.addcritpath, cb=cb)
        solver.debug && println("  added $(ncuts) cuts.")

    end
    addlazycallback(master.model, cbgbd)

    # solve master problem (has callback)
    status = solve(master.model, suppress_warnings=true)

    # TODO: extract bounds, nnodes?
    primal = Inf  # could compute from bestsol?
    dual = getobjectivevalue(master.model)  # correct bound?
    nnodes = -1
    Result(status, bestsol, primal, dual, nnodes)
end
