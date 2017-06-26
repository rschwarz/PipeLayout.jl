import PipeLayout: pwl_ineqs, reorient_fwdflow, pwl_inverse, pipesplit

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
    C = inst.ploss_coeff .* q .* abs.(q)

    model = Model(solver=solver.solver)

    # node positions
    @variable(model, xmin ≤ x[1:nnodes] ≤ xmax)
    @variable(model, ymin ≤ y[1:nnodes] ≤ ymax)

    # differences
    @variable(model, xmin - xmax ≤ dx[1:narcs] ≤ xmax - xmin)
    @variable(model, ymin - ymax ≤ dy[1:narcs] ≤ ymax - ymin)

    # squared pressure on nodes
    @variable(model, πlb[v] ≤ π[v=1:nnodes] ≤ πub[v])

    # auxiliary variable for pipe cost
    @variable(model, t[1:narcs] ≥ 0)

    # auxiliary variable for second-order cones
    @variable(model, L[1:narcs] ≥ 0)

    # fix terminal positions
    @constraint(model, fix[t=1:nterms], x[termidx[t]] == terms[t].x)
    @constraint(model, fiy[t=1:nterms], y[termidx[t]] == terms[t].y)

    # compute differences
    @constraint(model, diffx[a=1:narcs], dx[a] == x[T[a]] - x[H[a]])
    @constraint(model, diffy[a=1:narcs], dy[a] == y[T[a]] - y[H[a]])

    # second order cone for length
    @constraint(model, soc[a=1:narcs], dx[a]^2 + dy[a]^2 ≤ L[a]^2)

    # inequalities for PWL function
    ie = pwl_ineqs(Dm5, c)
    @assert size(ie) == (ndiams - 1, 3)
    @constraint(model, pwl[a=1:narcs, s=1:ndiams-1],
                ie[s,1]/C[a]*(π[T[a]] - π[H[a]]) + ie[s,2]*t[a] ≥ ie[s,3]*L[a])

    # lower bound for diameter impact (y)
    Dmax = inst.diameters[end].value
    @constraint(model, ylb[a=1:narcs], Dmax^5/C[a]*(π[T[a]] - π[H[a]]) ≥ L[a])

    # lower bound for cost factor (instead of upper bound on y)
    cmin = inst.diameters[1].cost
    @constraint(model, yub[a=1:narcs], t[a] ≥ cmin * L[a])

    # minimize total pipe cost
    @objective(model, :Min, sum(t[a] for a=1:narcs))

    model, x, y, t, π
end

function PipeLayout.optimize(inst::Instance, topo::Topology, solver::SOC)
    model, x, y, t, π = make_soc(inst, topo, solver)
    status = solve(model, suppress_warnings=true)
    objval = getobjectivevalue(model)
    nodes = map(Node, zip(getvalue(x), getvalue(y)))

    # compute l from t
    tsol = getvalue(t)
    L = pipelengths(Topology(nodes, topo.arcs))
    D, c = [d.value for d in inst.diameters], [d.cost for d in inst.diameters]
    zsol = tsol ./ L
    pwl_xs, pwl_ys = reverse(D.^(-5)), reverse(c) # must be sorted
    ysol = [pwl_inverse(pwl_xs, pwl_ys, z) for z in zsol]
    lsol = vcat([pipesplit(inst.diameters, y^(-1/5))' for y in ysol]...)
    @assert size(lsol) == (length(L), length(D))

    sol = Solution(nodes, lsol, getvalue(π))
    Result(status, sol, objval)
end
