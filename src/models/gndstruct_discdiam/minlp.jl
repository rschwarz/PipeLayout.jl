struct MINLP <: GroundStructureSolver
    solver             # to solve MINLP model
    contz::Bool
    timelimit::Float64 # seconds
    writemodels::Bool

    function MINLP(solver; contz=false, timelimit=Inf, writemodels=false)
        new(solver, contz, timelimit, writemodels)
    end
end

"""
Build MINLP model for instance on given ground structure.

To be solved with given MPB solver (capable of NLP).
Use continuous variables 0 ≤ z ≤ 1 if `contz` is true.
"""
function make_minlp(inst::Instance, topo::Topology, optimizer; contz=false)
    nodes, nnodes = topo.nodes, length(topo.nodes)
    arcs, narcs = topo.arcs, length(topo.arcs)
    terms, nterms = inst.nodes, length(inst.nodes)
    termidx = termindex(nodes, terms)
    ndiams = length(inst.diameters)

    # demand for all nodes, including junctions
    dem = fill(0.0, nnodes)
    for t in 1:nterms
        dem[termidx[t]] = inst.demand[t]
    end

    # pressure bounds for all nodes
    πlb_min = minimum([b.lb for b in inst.pressure])^2
    πub_max = maximum([b.ub for b in inst.pressure])^2
    πlb = fill(πlb_min, nnodes)
    πub = fill(πub_max, nnodes)
    πlb[termidx] = [b.lb^2 for b in inst.pressure]
    πub[termidx] = [b.ub^2 for b in inst.pressure]

    # adjacency lists
    inarcs, outarcs = [Int[] for v in 1:nnodes], [Int[] for v in 1:nnodes]
    for a in 1:narcs
        tail, head = arcs[a]
        push!(inarcs[head], a)
        push!(outarcs[tail], a)
    end
    tail = [arcs[a].tail for a in 1:narcs]
    head = [arcs[a].head for a in 1:narcs]

    # "big-M" bound for flow on arcs
    qmax = 0.5 * sum(abs.(inst.demand))

    L = pipelengths(topo)
    c = [diam.cost for diam in inst.diameters]
    Dm5 = [diam.value^(-5) for diam in inst.diameters]
    C = inst.ploss_coeff * L

    # always use SCIP directly
    model = JuMP.direct_model(MOI.instantiate(optimizer))

    # select arcs from topology
    @variable(model, y[1:narcs], Bin)

    # select diameter on arcs
    if contz # continuous
        @variable(model, 0 ≤ z[1:narcs, 1:ndiams] ≤ 1)
    else # discrete
        @variable(model, z[1:narcs, 1:ndiams], Bin)
    end

    # flow through arcs
    @variable(model, 0 <= q[1:narcs] <= qmax)

    # squared pressure on nodes
    @variable(model, πlb[v] ≤ π[v=1:nnodes] ≤ πub[v])

    # mass flow balance at nodes
    @constraint(model, flowbalance[v=1:nnodes],
                sum(q[a] for a=inarcs[v]) - sum(q[a] for a=outarcs[v]) == dem[v])

    # allow flow only for active arcs
    @constraint(model, active[a=1:narcs], q[a] <= qmax*y[a])

    # (conditional) pressure loss constraint (awkward nonlinear formulation)
    @NLconstraint(model, ploss[a=1:narcs],
                  y[a]*(π[tail[a]] - π[head[a]]) ==
                  C[a]*sum(Dm5[i]*z[a,i] for i=1:ndiams)*q[a]^2 )

    # choose diameter on active arcs
    @constraint(model, choice[a=1:narcs], sum(z[a,i] for i=1:ndiams) == y[a])

    # minimize total construction cost
    @objective(model, Min, sum(c[i] * L[a] * z[a,i] for a=1:narcs for i=1:ndiams))

    return model, y, z, q, π
end

function PipeLayout.optimize(inst::Instance, topo::Topology, solver::MINLP)
    if solver.timelimit < Inf
        MOI.set(solver.solver, MOI.TimeLimitSec, solver.timelimit)
    end
    model, y, z, q, π = make_minlp(inst, topo, solver.solver, contz=solver.contz)
    JuMP.optimize!(model)
    status = JuMP.termination_status(model)

    if solver.writemodels
        MOI.writeproblem(internalmodel(model), "minlp.orig.cip")
    end

    if status == MOI.INFEASIBLE
        return Result(status, nothing, Inf, Inf, 0)
    end

    zsol = round.(Bool, JuMP.value.(z) .>= 0.5)
    qsol = JuMP.value.(q)
    bestsol = CandSol(zsol, qsol, qsol.^2)

    primal = JuMP.objective_value(model)
    dual = JuMP.objective_bound(model)

    Result(status, bestsol, primal, dual, 0)
end
