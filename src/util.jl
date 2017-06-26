using JuMP
using MathProgBase

export ɛ, settimelimit!, stilltime

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
        internal = internalmodel(model)
        if isa(internal, MathProgBase.AbstractMathProgModel)
            # model is already built, want to modify current params
            MathProgBase.setparameters!(internal, TimeLimit=limit)
        else
            # model not yet build, want to modify future params
            MathProgBase.setparameters!(solver, TimeLimit=limit)
        end
    end
end

"is there still enough until final limit?"
function stilltime(finaltime; buffer=timebuffer)
    time() + buffer < finaltime
end
