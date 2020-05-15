struct NLP <: JunctionLocationSolver
    optimizer          # for underlying NLP model
    timelimit::Float64 # seconds

    function NLP(solver; timelimit=Inf)
        new(solver, timelimit)
    end
end

"""
Direct implementation of nonconvex NLP model, using split-pipe cost.

Optionally force discrete diameter choice for every segment.
"""
function make_nlp(inst::Instance, topo::Topology, solver::NLP;
                  discdiam=false)
    nodes, nnodes = topo.nodes, length(topo.nodes)
    arcs, narcs = topo.arcs, length(topo.arcs)
    terms, nterms = inst.nodes, length(inst.nodes)
    termidx = termindex(nodes, terms)
    ndiams = length(inst.diameters)
    xmin, xmax = extrema(node.x for node in terms)
    ymin, ymax = extrema(node.y for node in terms)

    # pressure bounds for all nodes
    πlb_min = minimum([b.lb for b in inst.pressure])^2
    πub_max = maximum([b.ub for b in inst.pressure])^2
    πlb = fill(πlb_min, nnodes)
    πub = fill(πub_max, nnodes)
    πlb[termidx] = [b.lb^2 for b in inst.pressure]
    πub[termidx] = [b.ub^2 for b in inst.pressure]

    q = uniq_flow(inst, topo)

    c = [diam.cost for diam in inst.diameters]
    Dm5 = [diam.value^(-5) for diam in inst.diameters]
    C = inst.ploss_coeff .* q .* abs.(q)

    model = JuMP.Model(solver.optimizer)

    # node positions
    @variable(model, xmin ≤ x[1:nnodes] ≤ xmax)
    @variable(model, ymin ≤ y[1:nnodes] ≤ ymax)

    # length of arcs
    @variable(model, L[1:narcs] ≥ 0)

    # relative length of segments
    if discdiam
        @variable(model, l[1:narcs,1:ndiams], Bin)
    else
        @variable(model, l[1:narcs,1:ndiams] ≥ 0)
    end

    # squared pressure on nodes
    @variable(model, πlb[v] ≤ π[v=1:nnodes] ≤ πub[v])

    # objective proxy variable
    @variable(model, z ≥ 0)


    # fix terminal positions
    @constraint(model, fix[t=1:nterms], x[termidx[t]] == terms[t].x)
    @constraint(model, fiy[t=1:nterms], y[termidx[t]] == terms[t].y)

    # couple positions and length
    @constraint(model, xylength[a=1:narcs],
                L[a]^2 == (x[arcs[a].head] - x[arcs[a].tail])^2 + (y[arcs[a].head] - y[arcs[a].tail])^2)

    # convex combination of diameters
    @constraint(model, convexcomb[a=1:narcs], sum(l[a,i] for i=1:ndiams) == 1.0)

    # pressure loss
    @constraint(model, ploss[a=1:narcs], π[arcs[a].tail] - π[arcs[a].head] ==
                C[a] * L[a] * sum(Dm5[i]*l[a,i] for i=1:ndiams))

    # minimize total construction cost
    @constraint(model, objfun,
                z ≥ sum(c[i] * L[a] * l[a,i] for a=1:narcs for i=1:ndiams))
    @objective(model, Min, z)

    model, x, y, L, l, π
end

function PipeLayout.optimize(inst::Instance, topo::Topology, solver::NLP)
    model, x, y, L, l, π = make_nlp(inst, topo, solver)
    JuMP.optimize!(model)
    status = JuMP.termination_status(model)
    if status in [MOI.INFEASIBLE, MOI.INFEASIBLE_OR_UNBOUNDED]
        return Result(status, Solution([], zeros(0,0), []), Inf)
    end
    objval = JuMP.objective_value(model)
    nodes = map(Node, zip(JuMP.value.(x), JuMP.value.(y)))
    sol = Solution(nodes, JuMP.value.(l), JuMP.value.(π))
    Result(status, sol, objval)
end
