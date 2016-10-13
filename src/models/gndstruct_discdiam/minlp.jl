immutable MINLP <: GroundStructureSolver
    solver # to solve MINLP model
    debug::Bool
    contz::Bool
    timelimit::Float64 # seconds

    function MINLP(solver; debug=false, contz=false, timelimit=Inf)
        new(solver, debug, contz, timelimit)
    end
end

"""
Build MINLP model for instance on given ground structure.

To be solved with given MPB solver (capable of NLP).
Use continuous variables 0 ≤ z ≤ 1 if `contz` is true.
"""
function make_minlp(inst::Instance, topo::Topology, solver; contz=false)
    nodes, nnodes = topo.nodes, length(topo.nodes)
    arcs, narcs = topo.arcs, length(topo.arcs)
    terms, nterms = inst.nodes, length(inst.nodes)
    termidx = [findfirst(nodes, t) for t in terms]
    all(termidx .> 0) || throw(ArgumentError("Terminals not part of topology"))
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
    const qmax = 0.5 * sum(abs(inst.demand))

    L = pipelengths(topo)
    c = [diam.cost for diam in inst.diameters]
    Dm5 = [diam.value^(-5) for diam in inst.diameters]
    C = inst.ploss_coeff * L

    model = Model(solver=solver)

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
                sum{q[a], a=inarcs[v]} - sum{q[a], a=outarcs[v]} == dem[v])

    # allow flow only for active arcs
    @constraint(model, active[a=1:narcs], q[a] <= qmax*y[a])

    # (conditional) pressure loss constraint (awkward nonlinear formulation)
    @NLconstraint(model, ploss[a=1:narcs],
                  y[a]*(π[tail[a]] - π[head[a]]) ==
                  C[a]*sum{Dm5[i]*z[a,i], i=1:ndiams}*q[a]^2 )

    # choose diameter on active arcs
    @constraint(model, choice[a=1:narcs], sum{z[a,i], i=1:ndiams} == y[a])

    # minimize total construction cost
    @objective(model, :Min, sum{c[i] * L[a] * z[a,i], a=1:narcs, i=1:ndiams})

    return model, y, z, q, π
end

function optimize(inst::Instance, topo::Topology, solver::MINLP)
    MathProgBase.setparameters!(solver.solver, TimeLimit=solver.timelimit)
    model, y, z, q, π = make_minlp(inst, topo, solver.solver, contz=solver.contz)
    status = solve(model)

    zsol = round(Bool, getvalue(z) .>= 0.5)
    qsol = getvalue(q)
    bestsol = CandSol(zsol, qsol, qsol.^2)

    dual = getobjbound(internalmodel(model))

    Result(status, bestsol, dual, 0)
end
