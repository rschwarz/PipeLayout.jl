using Unitful

export ploss_coeff, ploss_coeff_nice, isunitless

function nounit(value)
    try
        uconvert(Unitful.NoUnits, value)
        true
    catch
        false
    end
end

# Source for notation and formulas:
# "Evaluating Gas Network Capacities", T. Koch et al., SIAM MO21, 2015
#
# Some constants are taken from that book, or the related source code of
# Lamatto++, see: http://www.mso.math.fau.de/edom/projects/lamatto.html

# pipe roughness, from REKO data (4390 * 0.012mm, 926 * 0.006mm, 166 * 0.002mm)
const k = 0.012u"mm"

# pipe (mean) diameter for friction factor
const D_m = 1.0u"m"

"friction factor by Nikoradse, for turbulent flow"
function nikoradse(diameter, roughness)
    (log10(diameter/roughness) + 1.138)^(-2)
end

# friction factor
const λ = nikoradse(D_m, k)

# universal gas constant (https://en.wikipedia.org/wiki/Gas_constant)
const R = 8.3144598u"J/K/mol"

# molar mass of methane
const molmass = 16.04u"g/mol"

# specific gas constant
const R_s = R/molmass

# critical pressure (default in Lamatto++)
const p_c = 46.4512u"bar"

# critical temperature (default in Lamatto++)
const T_c = 192.033u"K"

# mean pressure
const p_m = 60u"bar"

# mean temperature (default in Lamatto++)
const T_m = 283.15u"K"

"compressibility factor by AGA8"
function aga8(pressure, temperature)
    # "reduced" pressure and temperature
    p_r = pressure / p_c
    T_r = temperature / T_c
    z = 1.0 + 0.257 * p_r - 0.533 * p_r/T_r
end

# z_m: mean compressibility factor
const z_m = aga8(p_m, T_m)

# Pressure Loss equation ("Weymouth"), when we care about length and diameter
#   p_o^2 - p_i^2 = L/D^5 C q|q|
const weymouth = (π/4)^(-2) * λ * R_s * z_m * T_m
@assert !nounit(weymouth)

# alternative coefficient with same order of magnitude for nice, round numbers
const weymouth_nice = 3000u"m^2/s^2"
@assert nounit(weymouth/weymouth_nice) # same unit

# Unit-less Weymouth coefficient: Add types for length, diameters and flow, then
# divide by the expected type to get just a float
const ploss_coeff_unit = 1u"(km*m^(-5)*(kg/s)^2) / (bar^2)"
const ploss_coeff = weymouth * ploss_coeff_unit
const ploss_coeff_nice = weymouth_nice * ploss_coeff_unit
@assert nounit(ploss_coeff)
@assert nounit(ploss_coeff_nice)
