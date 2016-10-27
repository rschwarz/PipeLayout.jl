immutable SOC <: JunctionLocationSolver
    solver             # for underlying NLP model
    timelimit::Float64 # seconds

    function SOC(solver; timelimit=Inf)
        new(solver, timelimit)
    end
end

"Convex reformulation of NLP model with SOC"
function make_soc(inst::Instance, topo::Topology, solver::SOC)
    # TODO: move this to another function that calls this one
    topo = reorient_fwdflow(inst, topo)

    nodes, nnodes = topo.nodes, length(topo.nodes)
    arcs, narcs = topo.arcs, length(topo.arcs)
    terms, nterms = inst.nodes, length(inst.nodes)
    termidx = [findfirst(nodes, t) for t in terms]
    all(termidx .> 0) || throw(ArgumentError("Terminals not part of topology"))
    ndiams = length(inst.diameters)
    xmin, xmax = extrema(node.x for node in terms)
    ymin, ymax = extrema(node.y for node in terms)

    # tail and head for arcs
    T = [a.tail for a in arcs]
    H = [a.head for a in arcs]

    # pressure bounds for all nodes
    πlb_min = minimum([b.lb for b in inst.pressure])^2
    πub_max = maximum([b.ub for b in inst.pressure])^2
    πlb = fill(πlb_min, nnodes)
    πub = fill(πub_max, nnodes)
    πlb[termidx] = [b.lb^2 for b in inst.pressure]
    πub[termidx] = [b.ub^2 for b in inst.pressure]

    # need correct orientation (= positive flow)
    q = uniq_flow(inst, topo)
    @assert all(q .> 0)

    c = [diam.cost for diam in inst.diameters]
    Dm5 = [diam.value^(-5) for diam in inst.diameters]
    C = inst.ploss_coeff .* q .* abs(q)

    model = Model(solver=solver.solver)

    # node positions
    @variable(model, xmin ≤ x[1:nnodes] ≤ xmax)
    @variable(model, ymin ≤ y[1:nnodes] ≤ ymax)

    # squared pressure on nodes
    @variable(model, πlb[v] ≤ π[v=1:nnodes] ≤ πub[v])

    # auxiliary variable for pipe cost
    @variable(model, t[1:narcs] ≥ 0)


    # fix terminal positions
    @constraint(model, fix[t=1:nterms], x[termidx[t]] == terms[t].x)
    @constraint(model, fiy[t=1:nterms], y[termidx[t]] == terms[t].y)

    # inequalities for PWL function
    ineqs = pwl_ineqs(Dm5, c)
    @constraint(model, pwl[a=1:narcs, s=1:ndiams-1],
                (x[T[a]] - x[H[a]])^2 + (y[T[a]] - y[H[a]])^2
                ≤ TODO*(π[T[a]] - π[H[a]]) + TODO*t[a])

    # lower bound for diameter impact (y)
    TODO

    # lower bound for cost factor (instead of upper bound on y)
    TODO

    # minimize total pipe cost
    @objective(model, sum{t[a], a=1:narcs})
end
