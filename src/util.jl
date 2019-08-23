using MathOptInterface
using JuMP
using SCIP

const MOI = MathOptInterface

export ɛ, settimelimit!, stilltime, add_cons

# tolerance for numerical comparison
const ɛ = 1e-6

"numerically (strictly) nonzero"
nonzero(x::Real; atol=ɛ) = abs(x) > atol

"numerically less than"
isapproxlt(x, y; atol=ɛ) = isapprox(x, y, atol=atol) || x < y

"numerically more than"
isapproxgt(x, y; atol=ɛ) = isapprox(x, y, atol=atol) || x > y

"check that value is within bounds"
function isapproxin(x, lb, ub; atol=ɛ)
    isapproxgt(x, lb, atol=atol) && isapproxlt(x, ub, atol=atol)
end

"check the special ordered set (type 1) property"
function isSOS1(xs; atol=ɛ)
    all([isapproxgt(x, 0.0, atol=atol) for x in xs]) || return false
    sum(xs .>= atol) <= 1
end

"check the special ordered set (type 2) property"
function isSOS2(xs; atol=ɛ)
    all([isapproxgt(x, 0.0, atol=atol) for x in xs]) || return false
    nonzeros = xs .>= atol
    sum(nonzeros) <= 1 || sum(nonzeros) == 2 &&
        (sum(nonzeros[1:end-1] .& nonzeros[2:end]) == 1)
end


const timebuffer = 10.0 # seconds

"update timelimit (in seconds) for internal solver"
function settimelimit!(model::JuMP.Model, solver, limit)
    limit = max(limit, timebuffer) # at least buffer
    if limit < Inf
        MOI.set(JuMP.backend(model), MOI.TimeLimitSec, limit)
    end
end

"is there still enough until final limit?"
function stilltime(finaltime; buffer=timebuffer)
    time() + buffer < finaltime
end

"""Extract solution values for given JuMP variables, for use in callbacks."""
function SCIP.sol_values(o::SCIP.Optimizer,
                         vars::AbstractArray{JuMP.VariableRef},
                         sol::Ptr{SCIP.SCIP_SOL}=C_NULL)
    return SCIP.sol_values(o, [JuMP.index(v) for v in vars], sol)
end

"""
Add JuMP constraint to MOI backend, converting the types.

This is useful to add the result of @build_constraint directly.
Try not to commit type piracy by using an idiosyncratic function name.

TODO: add new macro to combine this with @build_constraint!
"""
function add_cons(optimizer::MOI.AbstractOptimizer,
                  constraint::JuMP.AbstractConstraint)
    f = JuMP.moi_function(constraint.func)
    s = constraint.set
    return MOI.add_constraint(optimizer, f, s)
end
