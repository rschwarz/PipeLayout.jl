using SIUnits
using SIUnits.ShortUnits

# check consistency of units in pressure loss
const L = 1.0km
const D = 1.0m
const q = 1.0kg/s
const p_i = 80.0Bar
const C = L/D^5*weymouth
# pi^2 - po^2 = C*q|q|
const π_i = p_i^2
const π_o = π_i - C*q*abs(q)

# pressure is dropping (and units compatible)
@test π_o < π_i

ratio = (π_i/π_o)^0.5
@test typeof(ratio) == Float64
@test_approx_eq_eps ratio 1.0 1e-6
