using JuMP
using SCIP

mutable struct Problem
    model
    x
    y
end

function make_model(nonlin=false)::Problem
    m = Model(solver=SCIP.Optimizer(display_verblevel=3))

    @variable(m, x[1:9], Bin)
    @variable(m, -1 ≤ y[1:9] ≤ 1)

    @constraint(m, c1[i=1:9], x[i] <= 3y[i])
    if nonlin
        @NLconstraint(m, sum(y[i] for i=1:9)^2 == prod(y[i] for i=1:9))
    end

    @objective(m, Max, sum(y[i] - i*x[i] for i=1:9))

    Problem(m, x, y)
end

function add_lazycons(prob::Problem)
    # some state for the lazy callback
    counter = 1

    function lazycb(cb)
        println("### lazycb is called! ($counter)")
        @lazyconstraint(cb, prob.x[counter] == prob.y[counter])

        counter = counter % 9 + 1
    end

    addlazycallback(prob.model, lazycb)
end

println("<<< First try without nonlinear constraints:")
prob = make_model(false)
add_lazycons(prob)
solve(prob.model)

println("<<< Now try with nonlinear constraints:")
prob = make_model(true)
add_lazycons(prob)
solve(prob.model)
