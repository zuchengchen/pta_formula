from numpy import pi

def STT(f, A, k):
    yr = 365.25 * 24.0 * 3600.0
    fyr = f * yr # f/(yr^{-1}) = f*yr
    coef = (1.0 + k**2.0)/(24.0 * pi**2.0 * f**3.0 * (1.0 + k**2 * fyr**(-2.0/3.0)))
    return coef * A**2.0 * fyr**(-4.0/3.0)


def SGeneral(f, A, k):
    yr = 365.25 * 24.0 * 3600.0
    fyr = f * yr # f/(yr^{-1}) = f*yr
    coef = (1.0 + k**2.0)/(24.0 * pi**2.0 * f**3.0 * (1.0 + k**2 * fyr**(-2.0/3.0)))
    return coef * A**2.0 * fyr**(-2.0)

def SInflation(f, Omega):
    H0 = 2.192711267238057e-18
    return 3.0 * H0**2/(16.0 * pi**4 * f**5) * Omega   