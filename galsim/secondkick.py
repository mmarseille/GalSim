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
This file implements the 'second-kick' profile used in the geometric approximation of
an atmospheric PSF.
"""

import numpy as np

import galsim

from . import _galsim
from .gsobject import GSObject


class SecondKick(GSObject):
    def __init__(self, lam, r0, L0=np.inf, kcrit=0, flux=1, scale_unit=galsim.arcsec,
                 do_delta=True, gsparams=None):
        if L0 > 1e10:
            L0 = 1e10
        if isinstance(scale_unit, str):
            scale_unit = galsim.angle.get_angle_unit(scale_unit)
        self._scale_unit = scale_unit
        scale = scale_unit/galsim.arcsec
        self._sbsk = _galsim.SBSecondKick(lam, r0, L0, kcrit, flux, scale, do_delta, gsparams)
        self._delta_amplitude = self._sbsk.getDeltaAmplitude()
        if do_delta:
            self._sbdelta = _galsim.SBDeltaFunction(self._delta_amplitude, gsparams=gsparams)
            # self._sbsk = _galsim.SBSecondKick(lam, r0, L0, kcrit, flux-self._delta_amplitude,
            #                                   scale, False, gsparams)
            self._sbsk = _galsim.SBSecondKick(lam, r0, L0, kcrit, flux-self._delta_amplitude,
                                              scale, do_delta, gsparams)
            self._sbp = _galsim.SBAdd([self._sbsk, self._sbdelta], gsparams=gsparams)
        else:
            self._sbp = self._sbsk

    @property
    def lam(self):
        return self._sbsk.getLam()

    @property
    def r0(self):
        return self._sbsk.getR0()

    @property
    def L0(self):
        return self._sbsk.getL0()

    @property
    def kcrit(self):
        return self._sbsk.getKCrit()

    @property
    def scale_unit(self):
        return self._scale_unit

    @property
    def do_delta(self):
        return self._sbsk.getDoDelta()

    @property
    def delta_amplitude(self):
        return self._sbsk.getDeltaAmplitude()

    def phasePower(self, kappa):
        return self._sbsk.phasePower(kappa)

    def structureFunction(self, rho):
        return self._sbsk.structureFunction(rho)

    def structureFunctionDirect(self, rho):
        return self._sbsk.structureFunctionDirect(rho)

_galsim.SBSecondKick.__getinitargs__ = lambda self: (
    self.getLam(), self.getR0(), self.getL0(), self.getKCrit(), self.getFlux(), self.getScale(),
    self.getDoDelta(), self.getGSParams())
_galsim.SBSecondKick.__getstate__ = lambda self: None
_galsim.SBSecondKick.__repr__ = lambda self: \
    "galsim._galsim.SBSecondKick(%r, %r, %r, %r, %r, %r, %r, %r)"%self.__getinitargs__()
