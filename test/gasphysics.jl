import SIUnits
using SIUnits.ShortUnits

import PipeLayout: weymouth

facts("unit consistency in pressure loss") do
    const L = 100.0km
    const D = 1.0m
    const q = 100.0kg/s
    const p_i = 80.0Bar
    const C = L/D^5*weymouth

    # pi^2 - po^2 = C*q|q|
    @fact p_i^2 - C*q*abs(q) --> less_than(p_i^2) "pressure should drop"
end
