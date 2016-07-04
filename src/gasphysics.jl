import SIUnits
using SIUnits.ShortUnits

export Bar, Barg, weymouth

# adding Bar and Barg (gauge) as units
const Bar = (10^5)Pa
const Barg = atm

# Source for notation and formulas:
# "Evaluating Gas Network Capacities", T. Koch et al., SIAM MO21, 2015

# pipe roughness, from REKO data (4390 * 0.012mm, 926 * 0.006mm, 166 * 0.002mm)
const k = 0.012mm

# pipe (mean) diameter for friction factor
const D_m = 1.0m

"friction factor by Nikoradse, for turbulent flow"
function nikoradse(diameter, roughness)
    (log10(diameter/roughness) + 1.138)^(-2)
end

# friction factor
const λ = nikoradse(D_m, k)

# universal gas constant (https://en.wikipedia.org/wiki/Gas_constant)
const R = 8.3144598J/K/mol

# molar mass of methane (TODO other mixture?)
const molmass = 16.04g/mol

# specific gas constant
const R_s = R/molmass

# critical pressure (default in Lamatto++)
const p_c = 46.4512Bar

# critical temperature (default in Lamatto++)
const T_c = 192.033K

# mean pressure (TODO find rationale)
const p_m = 60Bar

# mean temperature (default in Lamatto++)
const T_m = 283.15K

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
@assert !isa(weymouth, Real)

# Unit-less Weymouth coefficient: Add types for length, diameters and flow, then
# divide by the expected type to get just a float
const ploss_coeff = weymouth * (km*m^(-5)*(kg/s)^2) / (Bar^2)
@assert isa(ploss_coeff, Real)
