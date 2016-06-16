import SIUnits

# Source for notation and formulas:
# "Evaluating Gas Network Capacities", T. Koch et al., SIAM MO21, 2015

# Pressure Loss equation ("Weymouth")
#
# p_o^2 - p_i^2 = L Λ q|q|
#
# or:
#
# p_o^2 - p_i^2 = L/D^5 C q|q|
#
# with C = (π/4)^(-2) λ R_s z_m T_m
#
# and friction factor (Nikoradse) λ = (2 log10(D/k) + 1.138)^(-2)

# k: roughness, from REKO data (4390 * 0.012mm, 926 * 0.006mm, 166 * 0.002mm)
roughness = 0.012*Milli*Meter

# D: mean diameter for friction factor
diameter_mean = 1.0 Meter

# λ: using Nikoradse for turbulent flow
friction_factor = (log10(diameter_mean/roughness) + 1.138)^(-2)

# R_s: specific gas constant
spec_gas_constant = 0.0 # TODO

# z_m: mean compressibility factor
mean_compressibility = 0.0 # TODO

# T_m:
mean_temperature = 0.0 # TODO

# C: the "Weymouth coefficient", when we care for variable length and diameter
C = (π/4)^(-2) friction_factor * spec_gas_constant *
            mean_compressibility * mean_temperature
