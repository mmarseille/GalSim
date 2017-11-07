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
                 doDelta=False, gsparams=None):
        # We lose stability if L0 gets too large.  This should be close enough to infinity for
        # all practical purposes though.
        if L0 > 1e10:
            L0 = 1e10
        # Need _scale_unit for repr roundtriping.
        self._scale_unit = scale_unit
        scale = scale_unit/galsim.arcsec
        self._sbvk = _galsim.SBVonKarman(lam, r0, L0, flux, scale, doDelta, gsparams)
        self._deltaAmplitude = self._sbvk.getDeltaAmplitude()
        # Add in a delta function with appropriate amplitude if requested.
        if doDelta:
            self._sbdelta = _galsim.SBDeltaFunction(self._deltaAmplitude, gsparams=gsparams)
            # A bit wasteful maybe, but params should be cached so not too bad to recreate _sbvk?
            self._sbvk = _galsim.SBVonKarman(lam, r0, L0, flux-self._deltaAmplitude, scale,
                                             doDelta, gsparams)

            GSObject.__init__(self, _galsim.SBAdd([self._sbvk, self._sbdelta], gsparams=gsparams))
        else:
            GSObject.__init__(self, self._sbvk)

    @property
    def lam(self):
        return self._sbvk.getLam()

    @property
    def r0(self):
        return self._sbvk.getR0()

    @property
    def L0(self):
        return self._sbvk.getL0()

    @property
    def scale_unit(self):
        return self._scale_unit
        # Type conversion makes the following not repr-roundtrip-able, so we store init input as a
        # hidden attribute.
        # return galsim.AngleUnit(self._sbvk.getScale())

    @property
    def doDelta(self):
        return self._sbvk.getDoDelta()

    @property
    def deltaAmplitude(self):
        return self._deltaAmplitude

    @property
    def halfLightRadius(self):
        return self._sbvk.getHalfLightRadius()

    def structureFunction(self, rho):
        return self._sbvk.structureFunction(rho)

    def __eq__(self, other):
        return (isinstance(other, galsim.VonKarman) and
        self.lam == other.lam and
        self.r0 == other.r0 and
        self.L0 == other.L0 and
        self.flux == other.flux and
        self.scale_unit == other.scale_unit and
        self.doDelta == other.doDelta and
        self.gsparams == other.gsparams)

    def __hash__(self):
        return hash(("galsim.VonKarman", self.lam, self.r0, self.L0, self.flux, self.scale_unit,
                     self.doDelta, self.gsparams))

    def __repr__(self):
        out = "galsim.VonKarman(lam=%r, r0=%r, L0=%r"%(self.lam, self.r0, self.L0)
        out += ", flux=%r"%self.flux
        if self.scale_unit != galsim.arcsec:
            out += ", scale_unit=%r"%self.scale_unit
        if self.doDelta:
            out += ", doDelta=True"
        out += ", gsparams=%r"%self.gsparams
        out += ")"
        return out

    def __str__(self):
        return "galsim.VonKarman(lam=%r, r0=%r, L0=%r)"%(self.lam, self.r0, self.L0)

_galsim.SBVonKarman.__getinitargs__ = lambda self: (
    self.getLam(), self.getR0(), self.getL0(), self.getFlux(), self.getScale(),
    self.getDoDelta(), self.getGSParams())
_galsim.SBVonKarman.__getstate__ = lambda self: None
_galsim.SBVonKarman.__repr__ = lambda self: \
    "galsim._galsim.SBVonKarman(%r, %r, %r, %r, %r, %r, %r)"%self.__getinitargs__()
