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
"""@file base.py
This file implements the von Karman atmospheric PSF profile.  A version which has the underlying
turbulence power spectrum truncated below a given scale is also available as a correction when using
geometric shooting through an atmospheric PhaseScreenPSF.
"""

import numpy as np

import galsim

from . import _galsim
from .gsobject import GSObject


class VonKarman(GSObject):
    def __init__(self, lam, r0, L0=np.inf, flux=1, scale_unit=galsim.arcsec,
                 gsparams=None):
        # We lose stability if L0 gets too large.  This should be close enough to infinity for
        # all practical purposes though.
        if L0 > 1e10:
            L0 = 1e10
        scale = scale_unit/galsim.radians
        GSObject.__init__(
            self,
            _galsim.SBVonKarman(
                lam*1e-9,  # nm -> m
                r0,
                L0,
                flux,
                scale,
                gsparams
            )
        )

    @property
    def lam(self):
        return self.SBProfile.getLam()*1e9  # m -> nm

    @property
    def r0(self):
        return self.SBProfile.getR0()

    @property
    def L0(self):
        return self.SBProfile.getL0()

    @property
    def scale_unit(self):
        return galsim.AngleUnit(self.SBProfile.getScale())

    def structureFunction(self, rho):
        return self.SBProfile.structureFunction(rho)
