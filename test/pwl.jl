import PipeLayout: linecoefs, pwl_ineqs

facts("linear interpolation of 2d points") do
    # simple line through origin
    s, o = linecoefs(0, 0, 1, 2)
    @fact s --> roughly(2)
    @fact o --> roughly(0)

    # order of points doesn't matter
    s, o = linecoefs(1, 2, 0, 0)
    @fact s --> roughly(2)
    @fact o --> roughly(0)

    # neg. slope with offset
    s, o = linecoefs(2, 2, 4, 1)
    @fact s --> roughly(-0.5)
    @fact o --> roughly(3)
end

facts("inequalities for convex PWL function") do
    # symmetric quadratic, with flat segment in the middle
    xs = [-5, -3, -1, 1, 3, 5]
    ineqs = pwl_ineqs(xs, xs.^2)
    @fact size(ineqs) --> (5, 3)
    @fact ineqs[1,:] ./ ineqs[end,:]   --> roughly([-1, 1, 1])
    @fact ineqs[2,:] ./ ineqs[end-1,:] --> roughly([-1, 1, 1])
    @fact ineqs[3,:] ./ ineqs[3,2]     --> roughly([ 0, 1, 1])
end
