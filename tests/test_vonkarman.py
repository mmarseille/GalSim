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
    # Check against pure python calculation
    vk = galsim.VonKarman(lam=500.0, r0=0.15, L0=25.0)
    np.testing.assert_approx_equal(vk.deltaAmplitude, 6.13892754335e-190)

    lam = 500.0
    r0 = 0.2
    L0 = 1e10
    vk = galsim.VonKarman(lam, r0, L0, scale_unit=galsim.arcsec)
    kolm = galsim.Kolmogorov(lam=lam, r0=r0)

    for k in np.linspace(0, vk.maxK(), 5):
        np.testing.assert_allclose(kolm.kValue(0,k).real, vk.kValue(0,k).real, rtol=0, atol=1e-5)


    # # Check which VonKarman profiles are actually constructible.
    # for lam in [300.0, 500.0, 1100.0]:
    #     for r0_500 in [0.05, 0.1, 0.2, 0.3  ]:
    #         r0 = r0_500*(lam/500)**(6./5)
    #         for L0 in [1e10, 100.0, 25.0, 10.0]:
    #             print("Attempting to use vk with lam, r0, L0 = {}".format((lam, r0, L0)))
    #             try:
    #                 vk = galsim.VonKarman(lam, r0, L0)
    #                 print(vk.stepK(), vk.maxK())
    #                 vk.drawImage()
    #                 vk.drawImage(method='phot', n_photons=100)
    #             except RuntimeError:
    #                 print("Failed")


if __name__ == "__main__":
    test_vk()
