# Copyright (c) 2012-2017 by the GalSim developers team on GitHub
# https://github.com/GalSim-developers
#
# This file is part of GalSim: The modular galaxy image simulation toolkit.
# https://github.com/GalSim-developers/GalSim
#
# GalSim is free software: redistribution and use in source and binary forms,
# with or without modification, are permitted provided that the following
# conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions, and the disclaimer given in the accompanying LICENSE
#    file.
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions, and the disclaimer given in the documentation
#    and/or other materials provided with the distribution.
#
"""@file dcr.py
Implementation of atmospheric differential chromatic refraction.

This file defines functions that return the refraction angle (the angle between the true and
apparent zenith angles of an object), as a function of zenith angle, wavelength, temperature,
pressure, and water vapor content.
"""

import galsim

def air_refractive_index_minus_one(wave, pressure=69.328, temperature=293.15, H2O_pressure=1.067):
    """Return the refractive index of air as function of wavelength.

    Uses the formulae given in Filippenko (1982), which appear to come from Edlen (1953),
    and Coleman, Bozman, and Meggers (1960).  The units of the original formula are non-SI,
    being mmHg for pressure (and water vapor pressure), and degrees C for temperature.  This
    function accepts SI units, however, and transforms them when plugging into the formula.

    The default values for temperature, pressure and water vapor pressure are expected to be
    appropriate for LSST at Cerro Pachon, Chile, but they are broadly reasonable for most
    observatories.

    @param wave           Wavelength array in nanometers
    @param pressure       Air pressure in kiloPascals.
    @param temperature    Temperature in Kelvins.
    @param H2O_pressure   Water vapor pressure in kiloPascals.

    @returns the refractive index minus 1.
    """
    P = pressure * 7.50061683 # kPa -> mmHg
    T = temperature - 273.15 # K -> C
    W = H2O_pressure * 7.50061683 # kPa -> mmHg

    sigma_squared = 1.0 / (wave * 1.e-3)**2.0 # inverse wavenumber squared in micron^-2
    n_minus_one = (64.328 + (29498.1 / (146.0 - sigma_squared))
                   + (255.4 / (41.0 - sigma_squared))) * 1.e-6
    n_minus_one *= P * (1.0 + (1.049 - 0.0157 * T) * 1.e-6 * P) / (720.883 * (1.0 + 0.003661 * T))
    n_minus_one -= (0.0624 - 0.000680 * sigma_squared)/(1.0 + 0.003661 * T) * W * 1.e-6
    return n_minus_one

def get_refraction(wave, zenith_angle, **kwargs):
    """Compute the angle of refraction for a photon entering the atmosphere.

    Photons refract when transitioning from space, where the refractive index n = 1.0 exactly, to
    air, where the refractive index is slightly greater than 1.0.  This function computes the
    change in zenith angle for a photon with a given wavelength.  Output is a positive number of
    radians, even though the apparent zenith angle technically decreases due to this effect.

    @param wave          Wavelength array in nanometers
    @param zenith_angle  as an Angle
    @param kwargs        Keyword arguments to pass to air_refractive_index() to override default
                         pressure, temperature, and/or H2O_pressure.

    @returns the absolute value of change in zenith angle in radians.
    """
    nm1 = air_refractive_index_minus_one(wave, **kwargs)
    # The following line is equivalent to:
    # n_squared = (nm1 + 1)**2
    # r0 = (n_squared - 1.0) / (2.0 * n_squared)
    r0 = nm1 * (nm1+2) / 2.0 / (nm1**2 + 2*nm1 + 1)
    return r0 * zenith_angle.tan()

def zenith_parallactic_angles(obj_coord, zenith_coord=None, HA=None, latitude=None):
    """Compute the zenith angle and parallactic angle of a celestial coordinate, given either
    the celestial coordinate of the zenith, or equivalently, the hour angle of the coordinate and
    the latitude of the observer.  This is useful for the function ChromaticAtmosphere() in the
    galsim.chromatic module.

    @param obj_coord     A CelestialCoord object for which the zenith and parallactic
                         angles will be computed.
    @param zenith_coord  A CelestialCoord indicating the coordinates of the zenith.
    @param HA            The hour angle (as an Angle) of the coordinate for which the
                         zenith and parallactic angles will be computed.
    @param latitude      The observer's latitude, as an Angle.

    @returns the tuple `(zenith_angle, parallactic_angle)`, each of which is an Angle.
    """
    if zenith_coord is None:
        if HA is None or latitude is None:
            raise ValueError("Need to provide either zenith_coord or (HA, latitude) to"
                             +"zenith_parallactic_angles()")
        zenith_coord = galsim.CelestialCoord(HA + obj_coord.ra, latitude)
    zenith_angle = obj_coord.distanceTo(zenith_coord)
    NCP = galsim.CelestialCoord(0.0*galsim.degrees, 90*galsim.degrees)
    parallactic_angle = obj_coord.angleBetween(zenith_coord, NCP)
    return zenith_angle, parallactic_angle
