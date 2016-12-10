import PipeLayout: linecoefs, pwl_ineqs, pwl_inverse

@testset "linear interpolation of 2d points" begin
    # simple line through origin
    s, o = linecoefs(0, 0, 1, 2)
    @test s ≈ 2
    @test o ≈ 0

    # order of points doesn't matter
    s, o = linecoefs(1, 2, 0, 0)
    @test s ≈ 2
    @test o ≈ 0

    # neg. slope with offset
    s, o = linecoefs(2, 2, 4, 1)
    @test s ≈ -0.5
    @test o ≈ 3
end

@testset "inequalities for convex PWL function" begin
    # symmetric quadratic, with flat segment in the middle
    xs = [-5, -3, -1, 1, 3, 5]
    ineqs = pwl_ineqs(xs, xs.^2)
    @test size(ineqs) == (5, 3)
    @test ineqs[1,:] ./ ineqs[end,:]   ≈ [-1, 1, 1]
    @test ineqs[2,:] ./ ineqs[end-1,:] ≈ [-1, 1, 1]
    @test ineqs[3,:] ./ ineqs[3,2]     ≈ [ 0, 1, 1]
end

@testset "inverses from monotoneous PWL functions" begin
    @testset "monotoneous increasing" begin
        xs = Float64[1, 2, 3, 4]
        ys = xs.^2
        @test_throws AssertionError pwl_inverse(xs, ys, 0.5)
        @test pwl_inverse(xs, ys, 1.0) ≈ 1.0
        @test pwl_inverse(xs, ys, 2.5) ≈ 1.5
        @test pwl_inverse(xs, ys, 16.0) ≈ 4.0
        @test_throws AssertionError pwl_inverse(xs, ys, 16.5)
    end

    @testset "monotoneous decreasing" begin
        xs = Float64[1, 2, 3, 4]
        ys = -xs.^2
        @test_throws AssertionError pwl_inverse(xs, ys, -0.5)
        @test pwl_inverse(xs, ys, -1.0) ≈ 1.0
        @test pwl_inverse(xs, ys, -2.5) ≈ 1.5
        @test pwl_inverse(xs, ys, -16.0) ≈ 4.0
        @test_throws AssertionError pwl_inverse(xs, ys, -16.5)
    end

    @testset "with constant pieces" begin
        xs = Float64[1, 2, 3]
        ys = Float64[1, 1, 5]
        @test_throws AssertionError pwl_inverse(xs, ys, 0.5)
        @test pwl_inverse(xs, ys, 1) ≈ 1.0
        @test pwl_inverse(xs, ys, 3) ≈ 2.5
        @test pwl_inverse(xs, ys, 5) ≈ 3.0
        @test_throws AssertionError pwl_inverse(xs, ys, 5.5)
    end
end
