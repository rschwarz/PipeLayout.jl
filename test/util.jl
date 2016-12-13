import PipeLayout: nonzero, isapproxlt, isapproxgt, isapproxin, isSOS1,
       isSOS2, select_subset

@testset "numerical comparisons" begin
    @test nonzero(0.0) == false
    @test nonzero(1e-15) == false
    @test nonzero(-1e-15) == false
    @test nonzero(1e-3) == true
    @test nonzero(-1e-3) == true

    @test isapproxlt(1.5, 2.0)
    @test isapproxlt(2.0 + 1e-10, 2.0)
    @test !isapproxlt(2.1, 2.0)

    @test isapproxgt(2.5, 2.0)
    @test isapproxgt(2.0 - 1e-10, 2.0)
    @test !isapproxgt(1.9, 2.0)

    @test !isapproxin(0.5, 1.0, 2.0)
    @test isapproxin(1.0 - 1e-9, 1.0, 2.0)
    @test isapproxin(1.5, 1.0, 2.0)
    @test isapproxin(2.0 + 1e-9, 1.0, 2.0)
    @test !isapproxin(2.5, 1.0, 2.0)

    @test !isSOS1([-1., 0.0])
    @test isSOS1( [0.0, 0.0])
    @test isSOS1( [0.5, 0.0])
    @test isSOS1( [0.0, 0.5])
    @test !isSOS1([0.3, 0.5])

    @test !isSOS2([-1., 0.0, 0.0])
    @test isSOS2( [0.0, 0.0, 0.0])
    @test isSOS2( [0.5, 0.0, 0.0])
    @test isSOS2( [0.0, 0.5, 0.3])
    @test !isSOS2([0.3, 0.0, 0.5])
    @test !isSOS2([1.0, 0.5, 0.3])
end

@testset "random subset sampling" begin
    base = collect(1:6)
    sample0 = select_subset(base, 0)
    sample3 = select_subset(base, 3)
    sample6 = select_subset(base, 6)

    @test length(sample0) == 0
    @test length(sample3) == 3
    @test unique(sample3) == sample3
    @test length(sample6) == 6
    @test unique(sample6) == sample6
    @test_throws AssertionError select_subset(base, 9)
end
