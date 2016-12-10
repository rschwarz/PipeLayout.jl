import PipeLayout: nonzero, select_subset

@testset "numerical comparisons" begin
    @test nonzero(0.0) == false
    @test nonzero(1e-15) == false
    @test nonzero(-1e-15) == false
    @test nonzero(1e-3) == true
    @test nonzero(-1e-3) == true
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
