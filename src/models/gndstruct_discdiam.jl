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
"""
immutable Master
    model::Model
    y::Vector{Variable}
    z::Array{Variable,2}
    q::Vector{Variable}
end

"Data type for candidate solutions from master problem."
immutable CandSol
    zsol::Array{Bool,2}
    qsol::Vector{Float64}
end

"Data type for dual solution of subproblem"
immutable SubDualSol
    μ::Array{Float64}
    λl::Array{Float64}
    λu::Array{Float64}
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

    model, y, z, q
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

"Construct no-good cut on z variables."
function nogood(master::Master, cand::CandSol)
    @assert size(master.z) == size(cand.zsol)
    nvars = length(master.z)
    zact = (cand.zsol .> 0.5)
    nact = sum(zact)
    coef = 2.0*zact - 1.0
    @constraint(master.model, sum{coef[i]*master.z[i], i=1:nvars} <= nact - 1)
    1
end

"Construct all Benders cuts from the solution of a subproblem."
function gndstruct_discdiam_cuts(inst::Instance, topo::Topology, master::Master,
                                 cand::CandSol, sub::SubDualSol)
    ncuts = 0
    ncuts += nogood(master, cand)
    ncuts
end

"Iteration based implementation of GBD."
function gndstruct_discdiam_algorithm(inst::Instance, topo::Topology)
    const maxiter = 100

    # initialize
    master = Master(gndstruct_discdiam_master(inst, topo)...)
    dual, status = 0.0, :Unknown

    for iter=1:maxiter
        println("Iter $(iter)")

        # resolve (relaxed) master problem, build candidate solution
        status = solve(master.model)
        if status == :Infeasible
            println("  relaxed master is infeasible :-(")
            break
        elseif status != :Optimal
            error("Unexpected status: $(:status)")
        end
        cand = CandSol(getvalue(master.z), getvalue(master.q))
        dual = getobjectivevalue(master.model)
        println("  dual bound: $(dual)")

        # solve subproblem (from scratch, no warmstart)
        submodel, π, Δl, Δu, ploss, plb, pub = gndstruct_discdiam_sub(inst, topo, cand)
        substatus = solve(submodel)
        @assert substatus == :Optimal "Slack model is always feasible"
        totalslack = getobjectivevalue(submodel)
        if totalslack ≈ 0.0
            println("  found feasible solution :-)")
            break
        end

        dualsol = SubDualSol(getdual(ploss), getdual(plb), getdual(pub))

        # generate cuts and add to master
        ncuts = gndstruct_discdiam_cuts(inst, topo, master, cand, dualsol)
        println("  added $(ncuts) cuts.")
    end
end
