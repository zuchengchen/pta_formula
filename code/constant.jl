"""
Gravitational constant.
"""
const G = 6.67408e-11


"""
Speed of light.
"""
const c = 2.997924580e8


""" 
1 day in seconds
"""
const day = 24.0 * 3600.0


"""
Year in seconds.
"""
const yr = 365.25day


"""
Euler's number
"""
const e = Base.MathConstants.e


"""
e in log10 base.
"""
const log10e = log10(e)


"""
10 in e base.
"""
const ln10 = log(10.0)


"""
Light year in meter.
"""
const ly = yr*c


"""
Parsec in meter.
"""
const pc = 3.085677581491367e16


"""
Kpc in meter.
"""
const Kpc = pc * 1.0e3


"""
Mpc in meter.
"""
const Mpc = pc * 1.0e6


"""
Gpc in meter.
"""
const Gpc = pc * 1.0e9


"""
G * Msun in natural unit.
"""
const GMsun = 1.327124400e20  # measured more precisely than Msun alone!


"""
Msun in natural unit.
"""
const Msun = GMsun / G


"""
Radius of sun in natural unit.
"""
const Rsun = GMsun / (c^2)


"""
"Period" of sun in natural unit.
"""
const Tsun = GMsun / (c^3)


"""
Planck constant.
"""
const h = 6.62607004e-34


"""
Frequency of 1/year.
"""
const fyr = 1.0 / yr


"""
Astronomical unit.
"""
const AU = 1.495978707e11


"""
The erg is a unit of energy equal to 10−7 joules.
"""
const erg = 1e-7


"""
for DM variation design matrix
"""
const DM_K = 2.41e-16