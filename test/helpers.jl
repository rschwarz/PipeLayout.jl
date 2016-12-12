"numerically less than"
approx_less(x, y) = x ≈ y || x < y

"numerically more than"
approx_greater(x, y) = x ≈ y || x > y

"check that value is within bounds"
approx_in(x, lb, ub) = approx_greater(x, lb) && approx_less(x, ub)

"check the special ordered set (type 1) property"
function is_SOS1(xs)
    ɛ = 1e-5
    all([approx_greater(x, 0.0) for x in xs]) || return false
    sum(xs .>= ɛ) <= 1
end

"check the special ordered set (type 2) property"
function is_SOS2(xs)
    ɛ = 1e-5
    all([approx_greater(x, 0.0) for x in xs]) || return false
    nonzeros = xs .>= ɛ
    sum(nonzeros) <= 1 || sum(nonzeros) == 2 &&
        (sum(nonzeros[1:end-1] & nonzeros[2:end]) == 1)
end

@testset "meta tests for helpers" begin
    @test approx_less(1.5, 2.0)
    @test approx_less(2.0 + 1e-10, 2.0)
    @test !approx_less(2.1, 2.0)

    @test approx_greater(2.5, 2.0)
    @test approx_greater(2.0 - 1e-10, 2.0)
    @test !approx_greater(1.9, 2.0)

    @test !approx_in(0.5, 1.0, 2.0)
    @test approx_in(1.0 - 1e-9, 1.0, 2.0)
    @test approx_in(1.5, 1.0, 2.0)
    @test approx_in(2.0 + 1e-9, 1.0, 2.0)
    @test !approx_in(2.5, 1.0, 2.0)

    @test !is_SOS1([-1., 0.0])
    @test is_SOS1( [0.0, 0.0])
    @test is_SOS1( [0.5, 0.0])
    @test is_SOS1( [0.0, 0.5])
    @test !is_SOS1([0.3, 0.5])

    @test !is_SOS2([-1., 0.0, 0.0])
    @test is_SOS2( [0.0, 0.0, 0.0])
    @test is_SOS2( [0.5, 0.0, 0.0])
    @test is_SOS2( [0.0, 0.5, 0.3])
    @test !is_SOS2([0.3, 0.0, 0.5])
    @test !is_SOS2([1.0, 0.5, 0.3])
end

#
# test data
#
"create a small instance (and topology) with a single pipe"
function single_pipe(;length=100.0, flow=200.0)
    nodes = [Node(0,0), Node(length,0)]
    demand = [-flow, flow]
    press = fill(Bounds(60.0, 80.0), size(nodes))
    diams = [Diameter(t...) for t in [(0.8, 1.0),(1.0, 1.2),(1.2, 1.5)]]
    inst = Instance(nodes, demand, press, diams)
    topo = Topology(nodes, [Arc(1,2)])
    inst, topo
end
