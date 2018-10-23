using Unitful

import PipeLayout: weymouth

@testset "unit consistency in pressure loss" begin
    L = 100.0u"km"
    D = 1.0u"m"
    q = 100.0u"kg/s"
    p_i = 80.0u"bar"
    C = L/D^5*weymouth

    # pi^2 - po^2 = C*q|q|,  pressure should drop
    @test p_i^2 - C*q*abs.(q) <= p_i^2
end
