import PipeLayout: optimize
using JuMP

export PipeDimLP, Solution, Result

immutable Solution
    lsol::Array{Float64,2}
    πsol::Vector{Float64}
end

immutable Result
    status::Symbol
    sol::Solution
    value::Float64
end

immutable PipeDimLP <: PipeDimensioningSolver
    lpsolver
end

function make_model(inst::Instance, topo::Topology, solver)
    nnodes = length(inst.nodes)
    length(topo.nodes) == nnodes ||
        throw(ArgumentError("Steiner nodes not allowed"))
    arcs, narcs = topo.arcs, length(topo.arcs)
    ndiams = length(inst.diameters)

    q = uniq_flow(inst, topo)

    model = Model(solver=solver)

    # squared pressure variables at nodes
    lb = [b.lb^2 for b in inst.pressure]
    ub = [b.ub^2 for b in inst.pressure]
    @variable(model, lb[v] <= π[v=1:nnodes] <= ub[v])

    # relative length of pipe segments with specific diameter
    @variable(model, 0.0 <= l[1:narcs, 1:ndiams] <= 1.0)

    # convex combination of pipe segments
    @constraint(model, convexcomb[a=1:narcs], sum(l[a,i] for i=1:ndiams) == 1.0)

    # pressure loss constraint, with variable factor from diameter choice
    L = pipelengths(topo)
    D5 = [diam.value^5 for diam in inst.diameters]
    C = inst.ploss_coeff * L .* q .* abs(q)
    @constraint(model, ploss[a=1:narcs], π[arcs[a].tail] - π[arcs[a].head]
                == sum(C[a] / D5[i] * l[a,i] for i=1:ndiams))

    # minimize total construction cost
    cost = [diam.cost for diam in inst.diameters]
    @objective(model, :Min, sum(cost[i] * L[a] * l[a,i] for a=1:narcs for i=1:ndiams))

    model, π, l
end

function optimize(inst::Instance, topo::Topology, solver::PipeDimLP)
    model, π, l = make_model(inst, topo, solver.lpsolver)
    status = solve(model)
    objval = status == :Optimal ? getobjectivevalue(model) : Inf
    sol = Solution(getvalue(l), getvalue(π))
    Result(status, sol, objval)
end
