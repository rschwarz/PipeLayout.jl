using JuMP

function pipedim_model(inst::Instance, topo::Topology)
    n = length(inst.nodes)
    length(topo.nodes) == n || throw(ArgumentError("Steiner nodes not allowed"))
    A, m = topo.arcs, length(topo.arcs)
    d = length(inst.diameters)

    q = uniq_flow(inst, topo)

    # TODO: think about scaling of variables and units
    model = Model()

    # squared pressure variables at nodes
    lb, ub = [b.lb for b in inst.pressure], [b.ub for b in inst.pressure]
    @variable(model, lb[v] <= π[v=1:n] <= ub[v])

    # relative length of pipe segments with specific diameter
    @variable(model, 0.0 <= l[1:m, 1:d] <= 1.0)

    # convex combination of pipe segments
    @constraint(model, convexcomb[a=1:m], sum{l[a,i], i=1:d} == 1.0)

    # pressure loss constraint, with variable factor from diameter choice
    L = pipelengths(topo)
    D = [diam.value^(-5) for diam in inst.diameters]
    # TODO: use weymouth constant and properly rm unit
    C = 1e-15 * L .* q .* abs(q)
    @constraint(model, ploss[a=1:m],
                π[A[a].tail] - π[A[a].head] == sum{C[a] * D[i] * l[a,i], i=1:d})

    # minimize total construction cost
    cost = [diam.cost for diam in inst.diameters]
    @objective(model, :Min, sum{cost[i] * L[a] * l[a,i], a=1:m, i=1:d})

    model, π, l
end
