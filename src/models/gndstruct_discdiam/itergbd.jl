"""
GBD iterations for Ground Structure with Discrete Diameters.

Solver object to store parameter values.
"""
struct IterGBD <: GroundStructureSolver
    addnogoods::Bool
    addcritpath::Bool
    maxiter::Int
    timelimit::Float64 # seconds
    writemodels::Bool
    mastersolver
    subsolver

    function IterGBD(mastersolver, subsolver;
                     addnogoods=false, addcritpath=true, maxiter::Int=100,
                     timelimit=Inf, writemodels=false)
        new(addnogoods, addcritpath, maxiter, timelimit, writemodels,
            mastersolver, subsolver)
    end
end

"Data type for dual solution of subproblem"
struct SubDualSol
    μ::Array{Float64}
    λl::Array{Float64}
    λu::Array{Float64}
end

"""
Data type for master problem.

Variable meanings:
  y: select arc
  z: select diameter for arc
  q: flow through arc
  ϕ: squared flow through arc (ϕ = q²)
"""
struct Master
    model::JuMP.Model
    y::Vector{JuMP.VariableRef}
    z::Array{JuMP.VariableRef, 2}
    q::Vector{JuMP.VariableRef}
    ϕ::Vector{JuMP.VariableRef}
end

"Build model for master problem (ground structure with discrete diameters)."
function make_master(inst::Instance, topo::Topology, optimizer)
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

    # adjacency lists
    inarcs, outarcs = [Int[] for v in 1:nnodes], [Int[] for v in 1:nnodes]
    for a in 1:narcs
        tail, head = arcs[a]
        push!(inarcs[head], a)
        push!(outarcs[tail], a)
    end

    # "big-M" bound for flow on arcs
    maxflow = 0.5 * sum(abs.(inst.demand))

    # always use direct mode for SCIP
    model = JuMP.direct_model(MOI.instantiate(optimizer))

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
                sum(q[a] for a=inarcs[v]) - sum(q[a] for a=outarcs[v]) == dem[v])

    # allow flow only for active arcs
    @constraint(model, active[a=1:narcs], q[a] <= maxflow*y[a])

    # choose diameter for active arcs
    @constraint(model, choice[a=1:narcs], sum(z[a,d] for d=1:ndiams) == y[a])

    # exclude antiparallel arcs.
    antiidx = antiparallelindex(topo)
    for a in 1:narcs
        if a < antiidx[a] # only at most once per pair
            @constraint(model, y[a] + y[antiidx[a]] ≤ 1)
        end
    end

    L = pipelengths(topo)
    c = [diam.cost for diam in inst.diameters]
    @objective(model, Min, sum(c[i] * L[a] * z[a,i] for a=1:narcs for i=1:ndiams))

    model, y, z, q, ϕ
end

"""
Build model for subproblem (ground structure with discrete diameters).

Corresponds to the domain relaxation with pressure loss overestimation, which
can be turned off via the flag `relaxed`.
"""
function make_sub(inst::Instance, topo::Topology, cand::CandSol, optimizer;
                  relaxed::Bool=true)
    nodes, nnodes = topo.nodes, length(topo.nodes)
    arcs, narcs = topo.arcs, length(topo.arcs)
    terms, nterms = inst.nodes, length(inst.nodes)
    termidx = termindex(nodes, terms)
    ndiams = length(inst.diameters)

    candarcs = filter(a -> any(cand.zsol[a,:]), 1:narcs)
    ncandarcs = length(candarcs)
    tail = [arcs[a].tail for a in candarcs]
    head = [arcs[a].head for a in candarcs]

    model = JuMP.Model(optimizer)

    # unconstrained variable for squared pressure, the bounds are added with
    # inequalities having slack vars.
    # for now, adding variables for all nodes, even if disconnected.
    @variable(model, π[1:nnodes])

    # overestimated pressure loss inequalities
    Dm5 = [diam.value^(-5) for diam in inst.diameters]
    C = inst.ploss_coeff * pipelengths(topo)[candarcs]
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

    @objective(model, Min, sum(Δl[v] + Δu[v] for v=1:nnodes))

    model, π, Δl, Δu, ploss, pres_lb, pres_ub
end

"Add tangent cut to quadratic inequality (ϕ ≥ q^2) if violated."
function quadratic_tangent(model, q, ϕ, qsol, ϕsol; cb=cb)
    violated = ϕsol < qsol^2 - ɛ
    if !violated
        return 0
    end

    # add 1st order Taylor approx:
    # q^2 ≈ 2*qsol*(q - qsol) + qsol^2 = 2*qsol*q - qsol^2
    if cb === nothing
        @constraint(model, ϕ ≥ 2*qsol*q - qsol^2)
    else
        @cb_constraint(cb, ϕ ≥ 2*qsol*q - qsol^2)
    end

    return 1
end

"""
Find tightest linear overestimator for given values v.

It's of the form: sum_i a_i x_i + sum_j b_j y_j + c
where x_i and y_j take values 0 or 1, with exactly one term active per sum.

The coefficients a, b, c are returned.
"""
function linear_overest(values::Matrix{Float64}, cand_i::Int, cand_j::Int,
                        optimizer)
    m, n = size(values)

    @assert 1 <= cand_i <= m
    @assert 1 <= cand_j <= n

    model = JuMP.Model(optimizer)

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
    @objective(model, Min, sum(t[i,j] for i=1:m for j=1:n))

    # solve it
    JuMP.optimize!(model)
    status = JuMP.termination_status(model)
    @assert status == MOI.OPTIMAL

    JuMP.value.(a), JuMP.value.(b), JuMP.value(c)
end

"Linearized & reformulated cut based on single path."
function pathcut(inst::Instance, topo::Topology, master::Master, cand::CandSol,
                 path::Vector{Arc}, solver; cb=nothing)
    aidx = arcindex(topo)
    pathidx = [aidx[arc] for arc in path]
    npath = length(path)
    ndiam = length(inst.diameters)
    zsol = cand.zsol[pathidx,:] # sparse solution
    qsol = cand.qsol[pathidx,:] # sparse solution
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
    ν = zsol * D.^(-5)
    @assert length(ν) == npath
    β = ν * (D.^5)'
    @assert size(β) == (npath, ndiam)
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
        β1st, β2nd = β[v-1, :], β[v, :]
        @assert size(β1st) == (ndiam,)
        @assert size(β2nd) == (ndiam,)
        βdiff = - repeat(β1st, 1, ndiam) + repeat(β2nd', ndiam, 1)
        @assert size(βdiff) == (ndiam, ndiam)
        supvalues[2:end, 2:end] = max.(βdiff * πub[node], βdiff * πlb[node])

        # the current values should be met exactly
        cand_i = findfirst(zsol[v-1,:])
        cand_j = findfirst(zsol[v,:])
        @assert cand_i ≠ 0 && cand_j ≠ 0

        # get coeffs of overestimation, assuming aux vars z_uv,0 and z_vw,0
        fix_i, fix_j = cand_i + 1, cand_j + 1
        cuv, cvw, c = linear_overest(supvalues, fix_i, fix_j, solver)

        # need to transform the coefficients to remove aux vars
        coeffs[v-1,:] += cuv[2:end] .- cuv[1]
        coeffs[v,:]   += cvw[2:end] .- cvw[1]
        offset += c + cuv[1] + cvw[1]
    end

    # - head of path
    head = path[end].head
    coeffs[end,:] += -1 * πlb[head] * β[end,:]

    # coefficients of ϕ
    C = inst.ploss_coeff * pipelengths(topo)[pathidx]
    α = ν .* C

    # add cut:  coeffs * z + offset ≥ α * ϕ
    z = master.z[pathidx,:]
    ϕ = master.ϕ[pathidx]

    if cb === nothing
        @constraint(
            master.model,
            sum(coeffs[a,i]*z[a,i] for a=1:npath for i=1:ndiam) + offset
            ≥ sum(α[a]*ϕ[a] for a=1:npath))
    else
        @cb_constraint(
            cb,
            sum(coeffs[a,i]*z[a,i] for a=1:npath for i=1:ndiam) + offset
            ≥ sum(α[a]*ϕ[a] for a=1:npath))
    end

    return 1
end

"Linearized & reformulated cuts based on critical paths."
function critpathcuts(inst::Instance, topo::Topology, master::Master,
                      cand::CandSol, sub::SubDualSol, solver; cb=nothing)
    ncuts = 0

    # compute dense dual flow
    narcs = length(topo.arcs)
    dualflow = fill(0.0, narcs)
    actarcs = vec(sum(cand.zsol, dims=2) .> 0)
    dualflow[actarcs] = sub.μ

    paths, pathflow = flow_path_decomp(topo, dualflow)

    # TODO: test that cuts are valid and separating
    for path in paths
        ncuts += pathcut(inst, topo, master, cand, path, solver, cb=cb)
    end

    for aidx in 1:narcs
        ncuts += quadratic_tangent(master.model, master.q[aidx], master.ϕ[aidx],
                                   cand.qsol[aidx], cand.ϕsol[aidx], cb=cb)
    end

    return ncuts
end

"Construct all Benders cuts from the solution of a subproblem."
function cuts(inst::Instance, topo::Topology, master::Master, cand::CandSol,
              sub::SubDualSol, solver; addnogoods=true, addcritpath=true)
    @assert any([addnogoods, addcritpath]) # must cut off!
    ncuts = 0
    if addnogoods
        ncuts += nogood(master.model, master.z, cand.zsol)
    end
    if addcritpath
        ncuts += critpathcuts(inst, topo, master, cand, sub, solver)
    end
    ncuts
end

"Iteration based implementation of GBD."
function PipeLayout.optimize(inst::Instance, topo::Topology, solver::IterGBD)
    run_gbd(inst, topo, solver.mastersolver, solver.subsolver,
            maxiter=solver.maxiter, timelimit=solver.timelimit,
            addnogoods=solver.addnogoods, addcritpath=solver.addcritpath,
            writemodels=solver.writemodels)
end

function run_gbd(inst::Instance, topo::Topology, mastersolver, subsolver;
                 maxiter::Int=100, timelimit=Inf, addnogoods=false,
                 addcritpath=true, writemodels=false)
    finaltime = time() + timelimit

    # initialize
    master = Master(make_master(inst, topo, mastersolver)...)
    dual, status = 0.0, MOI.OPTIMIZE_NOT_CALLED

    for iter=1:maxiter
        if !stilltime(finaltime)
            status = MOI.TIME_LIMIT
            @debug "Timelimit reached."
            return Result(status, nothing, Inf, dual, iter)
        end

        @debug "Iter $(iter)"

        # resolve (relaxed) master problem, build candidate solution
        writemodels && writeLP(master.model, "master_iter$(iter).lp", genericnames=false)
        settimelimit!(master.model, mastersolver, finaltime - time())
        JuMP.optimize!(master.model)
        status = JuMP.termination_status(master.model)
        if status == MOI.INFEASIBLE
            @debug "  relaxed master is infeasible :-("
            return Result(status, nothing, Inf, Inf, iter)
        elseif status in [MOI.TIME_LIMIT, MOI.NODE_LIMIT]
            return Result(status, nothing, Inf, dual, iter)
        elseif status == MOI.OPTIMAL
            # good, we continue below
        else
            error("Unexpected status: $(status)")
        end

        zsol = JuMP.value.(master.z)
        cand = CandSol(zsol .>= 0.5, JuMP.value.(master.q), JuMP.value.(master.ϕ))

        dual = JuMP.objective_value(master.model)
        @debug begin
            "  dual bound: $(dual)"
            "  cand. sol:$(Tuple.(findall(!iszero, cand.zsol)))"
        end

        # check whether candidate has tree topology
        ysol = JuMP.value.(master.y)
        candtopo = topology_from_candsol(topo, ysol)
        if !is_tree(candtopo)
            fullcandtopo = topology_from_candsol(topo, ysol, true)
            cycle = find_cycle(fullcandtopo)
            if length(cycle) == 0
                # TODO: Actually, this might be optimal, but it could also occur
                # when adding some irrelevant pipe is cheaper than increasing
                # the diameter. How to distinguish these cases?
                nogood(master.model, master.y, ysol)
                @debug "  skip disconnected topology with nogood."
            else
                avoid_topo_cut(master.model, master.y, topo, cycle)
                @debug "  skip non-tree topology, cycle: $(cycle)"
            end
            continue
        end

        # solve subproblem (from scratch, no warmstart)
        submodel, π, Δl, Δu, ploss, plb, pub = make_sub(inst, topo, cand, subsolver)
        writemodels && writeLP(submodel, "sub_relax_iter$(iter).lp", genericnames=false)
        settimelimit!(submodel, subsolver, finaltime - time())
        JuMP.optimize!(submodel)
        substatus = JuMP.termination_status(submodel)
        @assert substatus == MOI.OPTIMAL "Slack model is always feasible"
        totalslack = JuMP.objective_value(submodel)
        if totalslack ≈ 0.0
            # maybe only the relaxation is feasible, we have to check also the
            # "exact" subproblem with equations constraints.
            submodel2, _ = make_sub(inst, topo, cand, subsolver, relaxed=false)
            writemodels && writeLP(submodel2, "sub_exact_iter$(iter).lp", genericnames=false)
            settimelimit!(submodel2, subsolver, finaltime - time())
            JuMP.optimize!(submodel2)
            substatus2 = JuMP.termination_status(submodel2)
            @assert substatus2 == MOI.OPTIMAL "Slack model is always feasible"
            totalslack2 = JuMP.objective_value(submodel2)

            if totalslack2 ≈ 0.0
                @debug "  found feasible solution :-)"
                primal = JuMP.objective_value(master.model)
                return Result(substatus2, cand, primal, dual, iter)
            else
                # cut off candidate with no-good on z
                @debug "  subproblem/relaxation gap!"
                nogood(master.model, master.z, cand.zsol)
                continue
            end
        end

        dualsol = SubDualSol(JuMP.dual.(ploss), JuMP.dual.(plb), JuMP.dual.(pub))

        # generate cuts and add to master
        ncuts = cuts(inst, topo, master, cand, dualsol, subsolver,
                     addnogoods=addnogoods, addcritpath=addcritpath)
        @debug "  added $(ncuts) cuts."
    end

    Result(MOI.ITERATION_LIMIT, nothing, Inf, dual, maxiter)
end
