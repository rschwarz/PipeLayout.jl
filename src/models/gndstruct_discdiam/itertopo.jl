"""
Topology enumeration for Ground Structure with Discrete Diameters.

Solver object to store parameter values.
"""
immutable IterTopo <: GroundStructureSolver
    maxiter::Int
    timelimit::Float64 # seconds
    debug::Bool
    mastersolver
    subsolver
    writemodels::Bool

    function IterTopo(mastersolver, subsolver;
                      maxiter::Int=100, timelimit=Inf, debug=false,
                      writemodels=false)
        new(maxiter, timelimit, debug, mastersolver, subsolver, writemodels)
    end
end

"""
Build master model using only arc selection y and flows q.

This is used to enumerate (embedded) topologies, each of which is evaluated with
a "semi-subproblem". Returns model and variables y, q.
"""
function make_semimaster(inst::Instance, topo::Topology, solver)
    nodes, nnodes = topo.nodes, length(topo.nodes)
    arcs, narcs = topo.arcs, length(topo.arcs)
    terms, nterms = inst.nodes, length(inst.nodes)
    termidx = [findfirst(nodes, t) for t in terms]
    all(termidx .> 0) || throw(ArgumentError("Terminals not part of topology"))

    # demand for all nodes, including junctions
    dem = fill(0.0, nnodes)
    for t in 1:nterms
        dem[termidx[t]] = inst.demand[t]
    end

    # adjacency lists
    inarcs, outarcs = [Int[] for v in 1:nnodes], [Int[] for v in 1:nnodes]
    for a in 1:narcs
        tail, head = arcs[a]
        push!(inarcs[head], a)
        push!(outarcs[tail], a)
    end

    # "big-M" bound for flow on arcs
    const maxflow = 0.5 * sum(abs(inst.demand))

    model = Model(solver=solver)

    # select arcs from topology with y
    @variable(model, y[1:narcs], Bin)

    # flow through arcs
    @variable(model, 0 <= q[1:narcs] <= maxflow)

    # mass flow balance at nodes
    @constraint(model, flowbalance[v=1:nnodes],
                sum{q[a], a=inarcs[v]} - sum{q[a], a=outarcs[v]} == dem[v])

    # allow flow only for active arcs
    @constraint(model, active[a=1:narcs], q[a] <= maxflow*y[a])

    # exclude antiparallel arcs.
    antiidx = antiparallelindex(topo)
    for a in 1:narcs
        if a < antiidx[a] # only at most once per pair
            @constraint(model, y[a] + y[antiidx[a]] ≤ 1)
        end
    end

    # use cost of smallest diameter as topology-based relaxation
    L = pipelengths(topo)
    cmin = inst.diameters[1].cost
    @objective(model, :Min, sum{cmin * L[a] * y[a], a=1:narcs})

    model, y, q
end

"""
Build model for a problem where the topology y and flow q are fixed, but the
diameters z are still free. This is between the master- and subproblem and can
be used
 - as a primal heuristic
 - to generate stronger no-good cuts on the y variables

Returns model, list of candidate arcs and (sparse) variables z
"""
function make_semisub(inst::Instance, topo::Topology, cand::CandSol, solver)
    nnodes = length(topo.nodes)
    arcs = topo.arcs
    termidx = [findfirst(topo.nodes, t) for t in inst.nodes]
    ndiams = length(inst.diameters)
    candarcs = filter(a -> any(cand.zsol[a,:]), 1:length(arcs))
    ncandarcs = length(candarcs)
    tail = [arcs[a].tail for a in candarcs]
    head = [arcs[a].head for a in candarcs]

    πlb_min = minimum([b.lb for b in inst.pressure])^2
    πub_max = maximum([b.ub for b in inst.pressure])^2
    πlb = fill(πlb_min, nnodes)
    πub = fill(πub_max, nnodes)
    πlb[termidx] = [b.lb^2 for b in inst.pressure]
    πub[termidx] = [b.ub^2 for b in inst.pressure]

    L = pipelengths(topo)
    c = [diam.cost for diam in inst.diameters]
    Dm5 = [diam.value^(-5) for diam in inst.diameters]
    C = inst.ploss_coeff * L[candarcs] .* cand.qsol[candarcs].^2

    model = Model(solver=solver)

    @variable(model, z[1:ncandarcs, 1:ndiams], Bin)
    @variable(model, πlb[v] ≤ π[v=1:nnodes] ≤ πub[v])

    @constraint(model, ploss[a=1:ncandarcs],
                π[tail[a]] - π[head[a]] ≥ C[a] * sum{Dm5[i]*z[a,i], i=1:ndiams})
    @constraint(model, choice[a=1:ncandarcs], sum{z[a,i], i=1:ndiams} == 1)

    @objective(model, :Min, sum{c[i] * L[candarcs[a]] * z[a,i], a=1:ncandarcs, i=1:ndiams})

    model, candarcs, z
end

"Iteration based decomposition with semimaster and ~subproblem."
function optimize(inst::Instance, topo::Topology, solver::IterTopo)
    run_semi(inst, topo, solver.mastersolver, solver.subsolver,
             maxiter=solver.maxiter, timelimit=solver.timelimit,
             debug=solver.debug, writemodels=solver.writemodels)
end

function run_semi(inst::Instance, topo::Topology, mastersolver, subsolver;
                  maxiter::Int=100, timelimit=timelimit, debug=false,
                  writemodels=false)
    finaltime = time() + timelimit
    narcs = length(topo.arcs)
    ndiams = length(inst.diameters)

    # initialize
    primal, dual, status, bestsol = Inf, 0.0, :NotSolved, nothing
    mastermodel, y, q = make_semimaster(inst, topo, mastersolver)

    for iter=1:maxiter
        if !stilltime(finaltime)
            debug && println("Timelimit reached.")
            break
        end
        debug && println("Iter $(iter)")

        # resolve (relaxed) semimaster problem, build candidate solution
        writemodels && writeLP(mastermodel, "master_iter$(iter).lp")
        settimelimit!(mastermodel, mastersolver, finaltime - time())
        status = solve(mastermodel, suppress_warnings=true)
        if status == :Infeasible
            debug && println("  master problem is infeasible.")
            if primal == Inf
                # no solution was found
                return Result(:Infeasible, nothing, Inf, iter)
            else
                @assert bestsol ≠ nothing
                return Result(:Optimal, bestsol, dual, iter)
            end
        elseif status == :UserLimit
            return Result(:UserLimit, bestsol, dual, iter)
        elseif status == :Optimal
            # good, we continue below
        else
            error("Unexpected status: $(status)")
        end

        ysol, qsol = getvalue(y), getvalue(q)
        dual = getobjectivevalue(mastermodel)
        if debug
            println("  dual bound: $(dual)")
            println("  cand. sol:$(find(qsol))")
        end

        # stopping criterion
        if dual > primal - ɛ
            @assert bestsol ≠ nothing
            debug && println("  proved optimality of best solution.")
            return Result(:Optimal, bestsol, dual, iter)
        end

        # check whether candidate has tree topology
        candtopo = topology_from_candsol(topo, ysol)
        if !is_tree(candtopo)
            fullcandtopo = topology_from_candsol(topo, ysol, true)
            cycle = find_cycle(fullcandtopo)
            if length(cycle) == 0
                # TODO: see above (candidate disconnected?)
                nogood(mastermodel, y, ysol)
                debug && println("  skip disconnected topology with nogood.")
            else
                avoid_topo_cut(mastermodel, y, topo, cycle)
                debug && println("  skip non-tree topology, cycle: $(cycle)")
            end
            continue
        end

        # solve subproblem (from scratch, no warmstart)
        zcand = fill(false, narcs, ndiams)
        zcand[ysol .> 0.5, 1] = true
        cand = CandSol(zcand, qsol, qsol.^2)

        submodel, candarcs, z = make_semisub(inst, topo, cand, subsolver)
        writemodels && writeLP(submodel, "sub_iter$(iter).lp")
        settimelimit!(submodel, subsolver, finaltime - time())
        substatus = solve(submodel, suppress_warnings=true)
        if substatus == :Optimal
            # have found improving solution?
            newobj = getobjectivevalue(submodel)
            if newobj < primal
                primal = newobj
                znew = fill(false, narcs, ndiams)
                znew[candarcs,:] = (getvalue(z) .> 0.5)
                bestsol = CandSol(znew, qsol, qsol.^2)
                debug && println("  found improving solution: $(primal)")
            end
        elseif substatus != :Infeasible
            error("Unexpected status: $(:substatus)")
        end

        # generate nogood cut and add to master
        nogood(mastermodel, y, ysol)
    end

    Result(:UserLimit, bestsol, dual, maxiter)
end
