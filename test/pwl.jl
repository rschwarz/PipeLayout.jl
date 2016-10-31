import PipeLayout: linecoefs, pwl_ineqs, pwl_inverse

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

facts("inverses from monotoneous PWL functions") do
    context("monotoneous increasing") do
        xs = Float64[1, 2, 3, 4]
        ys = xs.^2
        @fact_throws pwl_inverse(xs, ys, 0.5)
        @fact pwl_inverse(xs, ys, 1.0) --> roughly(1.0)
        @fact pwl_inverse(xs, ys, 2.5) --> roughly(1.5) # between 1 and 4
        @fact pwl_inverse(xs, ys, 16.0) --> roughly(4.0)
        @fact_throws pwl_inverse(xs, ys, 16.5)
    end

    context("monotoneous decreasing") do
        xs = Float64[1, 2, 3, 4]
        ys = -xs.^2
        @fact_throws pwl_inverse(xs, ys, -0.5)
        @fact pwl_inverse(xs, ys, -1.0) --> roughly(1.0)
        @fact pwl_inverse(xs, ys, -2.5) --> roughly(1.5) # between -1 and -4
        @fact pwl_inverse(xs, ys, -16.0) --> roughly(4.0)
        @fact_throws pwl_inverse(xs, ys, -16.5)
    end

    context("with constant pieces") do
        xs = Float64[1, 2, 3]
        ys = Float64[1, 1, 5]
        @fact_throws pwl_inverse(xs, ys, 0.5)
        @fact pwl_inverse(xs, ys, 1) --> roughly(1.0)
        @fact pwl_inverse(xs, ys, 3) --> roughly(2.5)
        @fact pwl_inverse(xs, ys, 5) --> roughly(3.0)
        @fact_throws pwl_inverse(xs, ys, 5.5)
    end
end
