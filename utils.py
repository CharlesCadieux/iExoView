"""
Frequently used fonctions

@author: CharlesCadieux 2022
"""

import numpy as np


def find_nearest(array, value):
    """
    Return the index of the nearest value in array.
    """

    return (np.abs(array - value)).argmin()


def semiamp(P, Mp, Rp, Ms, i, e):
    """
    P: Orbital period (in days)
    Mp: Planetary mass (in Earth masses)
    Rp: Planetary radius (in Earth radii)
    Ms: Stellar mass (in solar masses)
    i: Orbital inclination (in degrees)
    e: Orbital eccentricity

    Return the radial velocity semi-amplitude K in m/s based on the input
    parameters.
    """

    G = 6.67408e-11 # Gravitational constant

    K = ((2 * np.pi * G / (P * 24 * 60 * 60)) ** (1 / 3) * (Mp * 5.972e24)
        * np.sin(np.radians(i)) / ((Ms * 1.9885e30)  + (Mp * 5.972e24))**(2 / 3)
        / np.sqrt(1 - e ** 2)) # in m/s

    return K


def mass_est(Rp):
    """
    Rp : Planetary radius (in Earth radii)

    Return the planetary mass estimate based on Chen and Kipping (2017) as
    implemented in Louie et al. (2018).
    """

    if Rp < 1.23:
        # Terrestrial size
        Mp = 0.9718 * Rp**3.58
        return Mp
    elif 1.23 <= Rp < 14.26:
        # Super-Earths to Jupiters size
        Mp = 1.436 * Rp**1.70
        return Mp
    else:
        return np.nan # Mass degeneracy for large radius


def radius_est(Mp):
    """
    Mp : Planetary mass (in Earth masses)

    Return the planetary radius estimate based on Chen and Kipping (2017).
    """

    if Mp < 2.04:
        # Terrestrial mass
        Rp = 1.008 * Mp**0.279
        return Rp
    elif 2.04 <= Mp < 131.6:
        # Neptunes mass
        Rp = 0.808 * Mp**0.589
        return Rp
    elif 131.6 <= Mp < 26637.2:
        # Jupiters mass
        Rp = 17.394 * Mp**(-0.04)
        return Rp
    else:
        # Stellar mass
        Rp = 1.476e-3 * Mp**0.88
        return Rp


def bulkrho(Mp, Rp):
    """
    Mp : Planetary mass (in Earth masses)
    Rp : Planet radius (in Earth radii)

    Return the mean density of the planet in g/cm3.
    """

    return (Mp * 5.972e24) / (4 * np.pi / 3 * (Rp * 6.371e6)**3) / 1000


def semima(P, Ms):
    """
    P: Orbital period (in days)
    Ms: Stellar mass (in solar masses)

    Return the semi-major axis of the planet orbit using Kepler's third law
    """

    return (Ms * (P / 365.25)**2)**(1/3)


def Teq(Rs, Teff, a):
    """
    Rs: Stellar radius (in solar radii)
    Teff: Stellar effective temperature (in K)
    a: Orbital separation (in au)

    Return the planetary equilibrium temperature for a null Bond albedo
    """

    return Teff * np.sqrt((Rs * 6.96340e8) / (a * 1.496e11)) * (0.25)**(0.25)


def atmosig(Teq, Rs, rho, mu):
    """
    Teq: Planetary equilibrium temperature (in K)
    Rs: Stellar radius (in solar radii)
    rho: Planet bulk density (in g/cm3)
    mu: Atmospheric mean molecular mass (in amu)

    Return the strength of transmission spectroscopy features in ppm for 5
    atmospheric scale height
    """

    G = 6.67408E-11 # Gravitational constant
    kb = 1.38064852E-23 # Boltzmann constant

    atmo = ((10 * kb * 3 / 4 / np.pi / G) * Teq / (Rs * 6.956E8)**2
        / (rho * 1000) / (mu * 1.66054E-27) * 1E6)

    return np.rint(atmo).astype(np.int32) # round to nearest integer (ppm)

def TSM(Rp, Mp, Rs, Teff, a, Jmag):
    """
    Rp : Planetary radius (in Earth radii)
    Mp : Planetary mass (in Earth masses)
    Rs: Stellar radius (in solar radii)
    Teff: Stellar effective temperature (in K)
    a: Orbital separation (in au)
    Jmag: Host star's J magnitude

    Return the Transmission Spectroscopy Metric from Kempton et al. 2018
    """

    if Rp < 1.5 : scale_fac = 0.190
    elif 1.5 <= Rp < 2.75 : scale_fac = 1.26
    elif 2.75 <= Rp < 4.0 : scale_fac = 1.28
    elif 4.0 <= Rp < 10.0 : scale_fac = 1.15
    else : scale_fac = 1.0

    teq = Teq(Rs, Teff, a)

    return scale_fac * Rp**3 * teq / Mp / Rs**2 * 10**(-Jmag/5)

def Planck(wav, T):
    """
    wav: Wavelength (in m)
    T: Temperature (in K)

    Return the blackbody spectral density (Planck law) in W/m2/sr/Hz
    """

    h = 6.62607004e-34 # Planck constant
    c = 2.99792458e+8 # speed of light
    kb = 1.38064852E-23 # Bolztmann constant
    a = 2.0 * h * c**2
    b = h * c / (wav * kb * T)

    return a / ((wav**5) * (np.exp(b) - 1.0))

def ESM(Rp, Rs, Teff, a, Kmag):
    """
    Rp : Planetary radius (in Earth radii)
    Rs: Stellar radius (in solar radii)
    Teff: Stellar effective temperature (in K)
    a: Orbital separation (in au)
    Kmag: Host star's K magnitude

    Return the Emission Spectroscopy Metric from Kempton et al. 2018
    """

    Tday = 1.10 * Teq(Rs, Teff, a) # Day side temperature estimate

    esm = (4.29e6 * Planck(7.5e-6, Tday) / Planck(7.5e-6, Teff)
        * ((Rp * 6.371e6) / (Rs * 6.96340e8))**2 * 10**(-Kmag/5))

    return esm
