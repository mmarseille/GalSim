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

import galsim

from . import _galsim
from .gsobject import GSObject


class SecondKick(GSObject):
    def __init__(self, lam, r0, L0, diam, obscuration, kcrit, flux=1, scale=1./206265, gsparams=None):
        self.lam = lam
        self.r0 = r0
        self.L0 = L0
        self.diam = diam
        self.obscuration = obscuration
        self.kcrit = kcrit
        self.scale = scale

        GSObject.__init__(
            self,
            _galsim.SBSecondKick(
                self.lam,
                self.r0,
                self.L0,
                self.diam,
                self.obscuration,
                self.kcrit,
                flux,
                self.scale,
                gsparams
            )
        )
        self._gsparams = gsparams
