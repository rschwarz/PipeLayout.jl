using SIUnits

# adding Bar and Barg (gauge) as units
const Bar = (1//10^5)Pascal
const Barg = 1Atmosphere

# Source for notation and formulas:
# "Evaluating Gas Network Capacities", T. Koch et al., SIAM MO21, 2015

# k: roughness, from REKO data (4390 * 0.012mm, 926 * 0.006mm, 166 * 0.002mm)
const roughness = 0.012*Milli*Meter

# D: mean diameter for friction factor
const diameter_mean = 1.0 Meter

# λ: using Nikoradse for turbulent flow
const friction_factor = (log10(diameter_mean/roughness) + 1.138)^(-2)

# R_s: specific gas constant (https://en.wikipedia.org/wiki/Gas_constant)
const spec_gas_constant = 8.3144598 Joule / Kelvin / Mole

# p_c: default critical pressure in Lamatto++
const crit_pressure = 46.4512 Bar

# T_c: default critical temperature in Lamatto++
const crit_temperature = 192.033 Kelvin

# z:
"compressibility factor by AGA8"
function compressibility(pressure, temperature)
    red_pressure = pressure / crit_pressure
    red_temperature = temperature / crit_temperature
    z = 1.0 + 0.257 red_pressure - 0.533 red_pressure / red_temperature
end

# p_m: TODO find rationale
const mean_pressure = 60 Bar

# T_m: default gas temperature in Lamatto++
const mean_temperature = 283.15 Kelvin

# z_m: mean compressibility factor
const mean_compressibility = compressibility(mean_pressure, mean_temperature)

# Pressure Loss equation ("Weymouth")
#
# p_o^2 - p_i^2 = L Λ q|q|
#
# or, because we care about length and diameter
#
# p_o^2 - p_i^2 = L/D^5 C q|q|
#
# with C = (π/4)^(-2) λ R_s z_m T_m
const C = ((π/4)^(-2) * friction_factor * spec_gas_constant
           * mean_compressibility * mean_temperature)
