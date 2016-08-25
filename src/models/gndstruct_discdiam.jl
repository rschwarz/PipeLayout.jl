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
function make_master(inst::Instance, topo::Topology;
                     solver=GLPKSolverMIP(msg_lev=0))
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

    # exclude antiparallel arcs.
    antiidx = antiparallelindex(topo)
    for a in 1:narcs
        if a < antiidx[a] # only at most once per pair
            @constraint(model, y[a] + y[antiidx[a]] ≤ 1)
        end
    end

    L = pipelengths(topo)
    c = [diam.cost for diam in inst.diameters]
    @objective(model, :Min, sum{c[i] * L[a] * z[a,i], a=1:narcs, i=1:ndiams})

    model, y, z, q, ϕ
end

"""
Build model for subproblem (ground structure with discrete diameters).

Corresponds to the domain relaxation with pressure loss overestimation, which
can be turned off via the flag `relaxed`.
"""
function make_sub(inst::Instance, topo::Topology, cand::CandSol;
                  solver=GLPKSolverLP(msg_lev=0), relaxed::Bool=true,
                  ploss_override::Float64=NaN)
    nodes, nnodes = topo.nodes, length(topo.nodes)
    arcs, narcs = topo.arcs, length(topo.arcs)
    terms, nterms = inst.nodes, length(inst.nodes)
    termidx = [findfirst(nodes, t) for t in terms]
    ndiams = length(inst.diameters)

    # TODO: why can't I set `ploss_coeff` as default value?
    if ploss_override === NaN # damn NaN
        ploss_override = ploss_coeff
    end

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
    C = ploss_override * pipelengths(topo)[candarcs]
    α = C .* cand.qsol[candarcs].^2 .* (cand.zsol[candarcs,:] * Dm5)
    if relaxed
        @constraint(model, ploss[ca=1:ncandarcs],
                    π[tail[ca]] - π[head[ca]] ≥ α[ca])
    else
        @constraint(model, ploss[ca=1:ncandarcs],
                    π[tail[ca]] - π[head[ca]] == α[ca])
    end

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

"""
Build master model using only arc selection y and flows q.

This is used to enumerate (embedded) topologies, each of which is evaluated with
a "semi-subproblem". Returns model and variables y, q.
"""
function make_semimaster(inst::Instance, topo::Topology;
                         solver=GLPKSolverMIP(msg_lev=0))
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
function make_semisub(inst::Instance, topo::Topology, cand::CandSol,
                      solver=GLPKSolverMIP(msg_lev=0))
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
    C = ploss_coeff * L[candarcs] .* cand.qsol[candarcs].^2

    model = Model(solver=solver)

    @variable(model, z[1:ncandarcs, 1:ndiams], Bin)
    @variable(model, πlb[v] ≤ π[v=1:nnodes] ≤ πub[v])

    @constraint(model, ploss[a=1:ncandarcs],
                π[tail[a]] - π[head[a]] ≥ C[a] * sum{Dm5[i]*z[a,i], i=1:ndiams})
    @constraint(model, choice[a=1:ncandarcs], sum{z[a,i], i=1:ndiams} == 1)

    @objective(model, :Min, sum{c[i] * L[a] * z[a,i], a=1:ncandarcs, i=1:ndiams})

    model, candarcs, z
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

"Add tangent cut to quadratic inequality (ϕ ≥ q^2) if violated."
function quadratic_tangent(model, q, ϕ, qsol, ϕsol)
    violated = ϕsol < qsol^2 - ɛ
    if !violated
        return 0
    end

    # add 1st order Taylor approx:
    # q^2 ≈ 2*qsol*(q - qsol) + qsol^2 = 2*qsol*q - qsol^2
    @constraint(model, ϕ ≥ 2*qsol*q - qsol^2)
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

    solver = GLPKSolverLP(msg_lev=0)
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
    @constraint(model, t[cand_i, cand_j] == 0)

    # be as tight as possible everywhere else
    @objective(model, :Min, sum{t[i,j], i=1:m, j=1:n})

    # solve it
    status = solve(model, suppress_warnings=true)
    @assert status == :Optimal

    getvalue(a), getvalue(b), getvalue(c)
end

"Linearized & reformulated cut based on single path."
function pathcut(inst::Instance, topo::Topology, master::Master, cand::CandSol,
                 path::Vector{Arc})
    aidx = arcindex(topo)
    pathidx = [aidx[arc] for arc in path]
    npath = length(path)
    ndiam = length(inst.diameters)
    zsol = cand.zsol[pathidx,:] # sparse solution
    D = [diam.value for diam in inst.diameters]

    # set loose pressure bounds for non-terminal nodes
    terminals = indexin(inst.nodes, topo.nodes)
    nnodes = length(topo.nodes)
    πlb_min = minimum([b.lb for b in inst.pressure])^2
    πub_max = maximum([b.ub for b in inst.pressure])^2
    πlb = fill(πlb_min, nnodes)
    πub = fill(πub_max, nnodes)
    πlb[terminals] = [b.lb^2 for b in inst.pressure]
    πub[terminals] = [b.ub^2 for b in inst.pressure]

    # coefficients of z in supremum expression
    β = (zsol * D.^(-5)) * (D.^5)'
    @assert all(β .> 0)

    # linearize the supremum, build up dense coefficient matrix
    coeffs = fill(0.0, (npath, ndiam))
    offset = 0.0
    @assert size(coeffs) == size(β)

    # - tail of path
    tail = path[1].tail
    coeffs[1,:] += πub[tail] * β[1,:]

    # - intermediate nodes, with in- and outgoing arcs
    for v in 2:npath
        uv, vw = path[v-1], path[v]
        @assert uv.head == vw.tail
        node = uv.head

        # prepare all possible values that need to be overestimated
        supvalues = zeros(ndiam + 1, ndiam + 1)
        # the case where no diameter is selected on the first arc
        supvalues[1, 2:end] = β[v, :] * πub[node]
        # similarly for no diameter on the second arc
        supvalues[2:end, 1] = - β[v-1, :] * πlb[node]
        # finally, when both arcs are active
        βdiff = repmat(β[v, :]', 1, ndiam) - repmat(β[v-1, :], ndiam, 1)
        supvalues[2:end, 2:end] = max(βdiff * πub[node], βdiff * πlb[node])

        # the current values should be met exactly
        cand_i = findfirst(zsol[v-1,:], true)
        cand_j = findfirst(zsol[v,:], true)
        @assert cand_i ≠ 0 && cand_j ≠ 0

        # get coeffs of overestimation, assuming aux vars z_uv,0 and z_vw,0
        cuv, cvw, c = linear_overest(supvalues, cand_i, cand_j)

        # need to transform the coefficients to remove aux vars
        coeffs[v-1,:] += cuv[2:end]' - cuv[1]
        coeffs[v,:]   += cvw[2:end]' - cvw[1]
        offset += c + cuv[1] + cvw[1]
    end

    # - head of path
    head = path[end].head
    coeffs[end,:] += πlb[head] * β[end,:]

    # add cut:  coeffs * z + offset ≥ α * ϕ
    z = master.z[pathidx,:]
    ϕ = master.ϕ[pathidx]
    C = ploss_coeff * pipelengths(topo)[pathidx]
    α = (zsol * D.^(-5)) .* C
    @constraint(master.model,
                sum{coeffs[a,i]*z[a,i], a=1:npath, i=1:ndiam} + offset ≥
                sum{α[a]*ϕ[a], a=1:npath})
    return 1
end

"Linearized & reformulated cuts based on critical paths."
function critpathcuts(inst::Instance, topo::Topology, master::Master,
                      cand::CandSol, sub::SubDualSol)
    ncuts = 0

    # compute dense dual flow
    narcs = length(topo.arcs)
    dualflow = fill(0.0, narcs)
    actarcs = collect(sum(cand.zsol, 2) .> 0)
    dualflow[actarcs] = sub.μ

    paths, pathflow = flow_path_decomp(topo, dualflow)

    # TODO: test that cuts are valid and separating
    for path in paths
        ncuts += pathcut(inst, topo, master, cand, path)
    end

    for aidx in 1:narcs
        ncuts += quadratic_tangent(master.model, master.q[aidx], master.ϕ[aidx],
                                   cand.qsol[aidx], cand.ϕsol[aidx])
    end

    return ncuts
end

"Cut off all y values for given, undirected topology"
function avoid_topo_cut(model, y, topo::Topology, edges::Vector{Arc})
    arcidx = arcindex(topo)
    antidx = antiparallelindex(topo)

    fwd = [arcidx[e] for e in edges]
    bwd = [antidx[a] for a in fwd]
    arcs = vcat(fwd, bwd)
    @constraint(model, sum{y[a], a in arcs} ≤ length(edges) - 1)
end

"Construct all Benders cuts from the solution of a subproblem."
function cuts(inst::Instance, topo::Topology, master::Master, cand::CandSol,
              sub::SubDualSol)
    # TODO: parametrize for computational experiments
    ncuts = 0
    ncuts += nogood(master.model, master.z, cand.zsol)
    ncuts += critpathcuts(inst, topo, master, cand, sub)
    ncuts
end

"Iteration based implementation of GBD."
function run(inst::Instance, topo::Topology; maxiter::Int=100, debug=false)

    # initialize
    master = Master(make_master(inst, topo)...)
    dual, status = 0.0, :NotSolved

    for iter=1:maxiter
        debug && println("Iter $(iter)")

        # resolve (relaxed) master problem, build candidate solution
        status = solve(master.model, suppress_warnings=true)
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
            cycle = find_cycle(candtopo)
            avoid_topo_cut(master.model, master.y, topo, cycle)
            debug && println("  skip non-tree topology, cycle: $(cycle)")
            continue
        end

        # solve subproblem (from scratch, no warmstart)
        submodel, π, Δl, Δu, ploss, plb, pub = make_sub(inst, topo, cand)
        substatus = solve(submodel, suppress_warnings=true)
        @assert substatus == :Optimal "Slack model is always feasible"
        totalslack = getobjectivevalue(submodel)
        if totalslack ≈ 0.0
            # TODO: should check whether subsolution is only feasible for
            # relaxation (overestimated pressure loss). If so, should cut off
            # current candidate with nogood. Maybe this is also a good time to
            # solve MIP-subproblem with free z?

            debug && println("  found feasible solution :-)")
            return Result(:Optimal, cand, dual, iter)
        end

        dualsol = SubDualSol(getdual(ploss), getdual(plb), getdual(pub))

        # generate cuts and add to master
        ncuts = cuts(inst, topo, master, cand, dualsol)
        debug && println("  added $(ncuts) cuts.")
    end

    Result(:UserLimit, nothing, dual, maxiter)
end

"Iteration based decomposition with semimaster and ~subproblem."
function run_semi(inst::Instance, topo::Topology; maxiter::Int=100, debug=false)
    narcs = length(topo.arcs)
    ndiams = length(inst.diameters)

    # initialize
    primal, dual, status, bestsol = Inf, 0.0, :NotSolved, nothing
    mastermodel, y, q = make_semimaster(inst, topo)

    for iter=1:maxiter
        debug && println("Iter $(iter)")

        # resolve (relaxed) semimaster problem, build candidate solution
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
        elseif status != :Optimal
            error("Unexpected status: $(:status)")
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
            debug && println("  skip non-tree topology")
            nogood(mastermodel, y, ysol)
            continue
        end

        # solve subproblem (from scratch, no warmstart)
        zcand = fill(false, narcs, ndiams)
        zcand[ysol .> 0.5, 1] = true
        cand = CandSol(zcand, qsol, qsol.^2)

        submodel, candarcs, z = make_semisub(inst, topo, cand)
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
