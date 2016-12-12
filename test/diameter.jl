import PipeLayout: serialmerge, pipesplit

@testset "split pipe utilities" begin
    diams = map(Diameter, [(1.0, 1.0), (1.2, 1.2), (1.4, 1.5)])

    @test serialmerge(diams, [1.0, 0.0, 0.0]) ≈ 1.0
    @test serialmerge(diams, [0.0, 1.0, 0.0]) ≈ 1.2
    @test serialmerge(diams, [0.0, 0.0, 1.0]) ≈ 1.4
    d = serialmerge(diams, [0.5, 0.5, 0.0])
    @test approx_in(d, 1.0, 1.2)

    @test pipesplit(diams, 1.0) ≈ [1.0, 0.0, 0.0]
    @test pipesplit(diams, 1.2) ≈ [0.0, 1.0, 0.0]
    @test pipesplit(diams, 1.4) ≈ [0.0, 0.0, 1.0]
    @test pipesplit(diams, d)   ≈ [0.5, 0.5, 0.0]
end
