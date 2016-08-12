import PipeLayout: digraph_from_topology
using GLPKMathProgInterface
using JuMP

export CandSol, SubDualSol, Master

"""
Data type for master problem.

Variable meanings:
  y: select arc
  z: select diameter for arc
  q: flow through arc
  ϕ: squared flow through arc (ϕ = q²)
"""
immutable Master
    model::Model
    y::Vector{Variable}
    z::Array{Variable,2}
    q::Vector{Variable}
    ϕ::Vector{Variable}
end

"Data type for candidate solutions from master problem."
immutable CandSol
    zsol::Array{Bool,2}
    qsol::Vector{Float64}
    ϕsol::Vector{Float64}
end

"Data type for dual solution of subproblem"
immutable SubDualSol
    μ::Array{Float64}
    λl::Array{Float64}
    λu::Array{Float64}
end

"Result data from GBD algorithm"
immutable Result
    status::Symbol
    solution
    dualbound::Float64
    niter::Int
end

function topology_from_candsol(topo::Topology, ysol::Vector{Float64})
    actarcs = topo.arcs[ysol .> 0.5]
    tails = [arc.tail for arc in actarcs]
    heads = [arc.head for arc in actarcs]
    nodeidx = tails ∪ heads
    nodemap = fill(0, length(topo.nodes))
    for (target, idx) in enumerate(nodeidx)
        nodemap[idx] = target
    end
    Topology(topo.nodes[nodeidx],
             [Arc(nodemap[arc.tail], nodemap[arc.head]) for arc in actarcs])
end

"Build model for master problem (ground structure with discrete diameters)."
function gndstruct_discdiam_master(inst::Instance, topo::Topology;
                                   solver=GLPKSolverMIP())
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

    # select single diameter for arc with z
    @variable(model, z[1:narcs, 1:ndiams], Bin)

    # flow through arcs
    @variable(model, 0 <= q[1:narcs] <= maxflow)
    @variable(model, 0 <= ϕ[1:narcs] <= maxflow^2)

    # secant cut for ϕ = q²
    @constraint(model, secant[a=1:narcs], ϕ[a] ≤ maxflow*q[a])

    # mass flow balance at nodes
    @constraint(model, flowbalance[v=1:nnodes],
                sum{q[a], a=inarcs[v]} - sum{q[a], a=outarcs[v]} == dem[v])

    # allow flow only for active arcs
    @constraint(model, active[a=1:narcs], q[a] <= maxflow*y[a])

    # choose diameter for active arcs
    @constraint(model, choice[a=1:narcs], sum{z[a,d], d=1:ndiams} == y[a])

    L = pipelengths(topo)
    c = [diam.cost for diam in inst.diameters]
    @objective(model, :Min, sum{c[i] * L[a] * z[a,i], a=1:narcs, i=1:ndiams})

    model, y, z, q, ϕ
end

"""
Build model for subproblem (ground structure with discrete diameters).

Corresponds to the domain relaxation with pressure loss overestimation.
"""
function gndstruct_discdiam_sub(inst::Instance, topo::Topology, cand::CandSol;
                                solver=GLPKSolverLP())
    nodes, nnodes = topo.nodes, length(topo.nodes)
    arcs, narcs = topo.arcs, length(topo.arcs)
    terms, nterms = inst.nodes, length(inst.nodes)
    termidx = [findfirst(nodes, t) for t in terms]
    ndiams = length(inst.diameters)

    candarcs = filter(a -> any(cand.zsol[a,:]), 1:narcs)
    ncandarcs = length(candarcs)
    tail = [arcs[a].tail for a in candarcs]
    head = [arcs[a].head for a in candarcs]

    model = Model(solver=solver)

    # unconstrained variable for squared pressure, the bounds are added with
    # inequalities having slack vars.
    # for now, adding variables for all nodes, even if disconnected.
    @variable(model, π[1:nnodes])

    # overestimated pressure loss inequalities
    Dm5 = [diam.value^(-5) for diam in inst.diameters]
    C = ploss_coeff * pipelengths(topo)[candarcs]
    α = C .* cand.qsol[candarcs].^2 .* (cand.zsol[candarcs,:] * Dm5)
    @constraint(model, ploss[ca=1:ncandarcs], π[tail[ca]] - π[head[ca]] ≥ α[ca])

    # slack variables to relax the lower and upper bounds for π.
    termlb = [b.lb^2 for b in inst.pressure]
    lb = fill(minimum(termlb), nnodes)
    lb[termidx] = termlb
    termub = [b.ub^2 for b in inst.pressure]
    ub = fill(maximum(termub), nnodes)
    ub[termidx] = termub

    @variable(model, Δl[1:nnodes] ≥ 0)
    @variable(model, Δu[1:nnodes] ≥ 0)

    @constraint(model, pres_lb[v=1:nnodes], π[v] + Δl[v] ≥ lb[v])
    @constraint(model, pres_ub[v=1:nnodes], π[v] - Δu[v] ≤ ub[v])

    @objective(model, :Min, sum{Δl[v] + Δu[v], v=1:nnodes})

    model, π, Δl, Δu, ploss, pres_lb, pres_ub
end

"Generic no-good cut"
function nogood{T<:Real}(model::Model, vars::Array{Variable}, sol::Array{T})
    @assert size(vars) == size(sol)
    nvars = length(vars)
    active = (sol .> 0.5)
    nactive = sum(active)
    coef = 2.0*active - 1.0
    @constraint(model, sum{coef[i]*vars[i], i=1:nvars} <= nactive - 1)
    return 1
end

"""
Find tightest linear overestimator for given values v.

It's of the form: sum_i a_i x_i + sum_j b_j y_j + c
where x_i and y_j take values 0 or 1, with exactly one term active per sum.

The coefficients a, b, c are returned.
"""
function linear_overest(values::Matrix{Float64}, cand_i::Int, cand_j::Int)
    m, n = size(values)

    @assert 1 <= cand_i <= m
    @assert 1 <= cand_j <= n

    solver = GLPKSolverLP()
    model = Model(solver=solver)

    # coefficients to be found
    @variable(model, a[1:m])
    @variable(model, b[1:n])
    @variable(model, c)

    # slack to minimize
    @variable(model, t[1:m,1:n] ≥ 0)

    # overestimate at all points
    @constraint(model, overest[i=1:m, j=1:n], a[i] + b[j] + c == values[i,j] + t[i,j])

    # be exact at candidate solution
    @constraint(model, t[cand_i, cand_j] ≤ 0)

    # be as tight as possible everywhere else
    @objective(model, :Min, sum{t[i,j], i=1:m, j=1:n})

    # solve it
    status = solve(model)
    @assert status == :Optimal

    getvalue(a), getvalue(b), getvalue(c)
end

"Linearized & reformulated cut based on single path."
function gndstruct_discdiam_pathcut(inst::Instance, topo::Topology,
                                    master::Master, cand::CandSol,
                                    path::Vector{Arc})
    aidx = arcindex(topo)
    pathidx = [aidx[arc] for arc in path]
    npath = length(path)
    ndiam = length(inst.diameters)
    zsol = cand.zsol[pathidx,:] # sparse solution
    D = [diam.value for diam in inst.diameters]
    πlb = [b.lb^2 for b in inst.bounds]
    πub = [b.ub^2 for b in inst.bounds]

    # coefficients of z in supremum expression
    β =  (zsol * D.^(-5)) * (D.^5)'
    @assert all(β .> 0)

    # linearize the supremum, build up dense coefficient matrix
    coeffs = fill(0.0, (npath, ndiam))
    @assert size(coeffs) == size(β)

    # - tail of path
    tail = path[1].tail
    coeffs[1,:] += πub[tail] * β[1,:]

    # - intermediate nodes
    # TODO: call linear_overest

    # - head of path
    head = path[end].head
    coeffs[end,:] += πlb[head] * β[end,:]

    # TODO: create cons  coeffs * z ≥ sum α ϕ
end

"Linearized & reformulated cuts based on critical paths."
function gndstruct_discdiam_critpathcuts(inst::Instance, topo::Topology,
                                         master::Master, cand::CandSol,
                                         sub::SubDualSol)
    # TODO: do path decomposition

    # TODO: add one cut for every path

    # TODO: add new linearizations for squared flows if separating
end

"Construct all Benders cuts from the solution of a subproblem."
function gndstruct_discdiam_cuts(inst::Instance, topo::Topology, master::Master,
                                 cand::CandSol, sub::SubDualSol)
    ncuts = 0
    ncuts += nogood(master.model, master.z, cand.zsol)
    ncuts
end

"Iteration based implementation of GBD."
function gndstruct_discdiam_algorithm(inst::Instance, topo::Topology;
                                      maxiter::Int=100, debug=false)

    # initialize
    master = Master(gndstruct_discdiam_master(inst, topo)...)
    dual, status = 0.0, :NotSolved

    for iter=1:maxiter
        debug && println("Iter $(iter)")

        # resolve (relaxed) master problem, build candidate solution
        status = solve(master.model)
        if status == :Infeasible
            debug && println("  relaxed master is infeasible :-(")
            return Result(:Infeasible, nothing, Inf, iter)
        elseif status != :Optimal
            error("Unexpected status: $(:status)")
        end
        cand = CandSol(getvalue(master.z), getvalue(master.q), getvalue(master.ϕ))
        dual = getobjectivevalue(master.model)
        if debug
            println("  dual bound: $(dual)")
            j,i,_ = findnz(cand.zsol')
            println("  cand. sol:$(collect(zip(i,j)))")
        end

        # check whether candidate has tree topology
        candtopo = topology_from_candsol(topo, getvalue(master.y))
        if !is_tree(candtopo)
            debug && println("  skip non-tree topology")
            nogood(master.model, master.y, getvalue(master.y))
            continue
        end

        # solve subproblem (from scratch, no warmstart)
        submodel, π, Δl, Δu, ploss, plb, pub = gndstruct_discdiam_sub(inst, topo, cand)
        substatus = solve(submodel)
        @assert substatus == :Optimal "Slack model is always feasible"
        totalslack = getobjectivevalue(submodel)
        if totalslack ≈ 0.0
            debug && println("  found feasible solution :-)")
            return Result(:Optimal, cand, dual, iter)
        end

        dualsol = SubDualSol(getdual(ploss), getdual(plb), getdual(pub))

        # generate cuts and add to master
        ncuts = gndstruct_discdiam_cuts(inst, topo, master, cand, dualsol)
        debug && println("  added $(ncuts) cuts.")
    end

    Result(:UserLimit, nothing, dual, maxiter)
end
