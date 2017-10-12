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

from __future__ import print_function
import numpy as np
import os
import sys

from galsim_test_helpers import *

try:
    import galsim
except ImportError:
    sys.path.append(os.path.abspath(os.path.join(path, "..")))
    import galsim


def test_vk():
    lam = 500.0
    r0 = 0.2
    L0 = 1e100
    kcrit = 0.0
    vk = galsim.VonKarman(lam, r0, L0, kcrit, scale_unit=galsim.arcsec)
    kolm = galsim.Kolmogorov(lam=lam, r0=r0)

    for k in [0.0, 1.0, 3.0, 10.0, 30.0]:
        np.testing.assert_allclose(kolm.kValue(0,k).real, vk.kValue(0,k).real, rtol=0, atol=3e-4)

if __name__ == "__main__":
    test_vk()
