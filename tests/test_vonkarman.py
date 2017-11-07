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


@timer
def test_vk():
    """Test the generation of VonKarman profiles
    """
    if __name__ == '__main__':
    # if __name__ != '__main__':
        lams = [300.0, 500.0, 1100.0]
        r0_500s = [0.05, 0.15, 0.3]
        L0s = [1e10, 25.0, 10.0]
        doDeltas = [False, True]
    else:
        lams = [500.0]
        r0_500s = [0.2]
        L0s = [25.0]
        doDeltas = [False]
    for lam in lams:
        for r0_500 in r0_500s:
            r0 = r0_500*(lam/500)**(6./5)
            for L0 in L0s:
                for doDelta in doDeltas:
                    kwargs = {'lam':lam, 'r0':r0, 'L0':L0, 'doDelta':doDelta}
                    print(kwargs)

                    vk = galsim.VonKarman(**kwargs)
                    np.testing.assert_almost_equal(vk.flux, 1.0)

                    check_basic(vk, "VonKarman")
                    do_pickle(vk)
                    do_pickle(vk.SBProfile)
                    do_pickle(vk.SBProfile, lambda x: (x.getFlux(), x.getGSParams()))

                    vk = galsim.VonKarman(**kwargs, flux=2.2)
                    np.testing.assert_almost_equal(vk.flux, 2.2)


@timer
def test_vk_delta():
    """Test a VonKarman with a significant delta-function amplitude"""
    kwargs = {'lam':1100.0, 'r0':0.8, 'L0':5.0, 'flux':2.2}
    vk = galsim.VonKarman(**kwargs)
    # This profile has more than 15% of its flux in the delta-function component.
    np.testing.assert_array_less(0.15, vk.deltaAmplitude/vk.flux)
    # If doDelta is False (the default), then the asymptotic kValue should still be zero.
    np.testing.assert_almost_equal(vk.kValue(1e10, 0).real, 0.0)
    # But if we use doDelta=True, then the asymptotic kValue should be that of the delta function.
    vkd = galsim.VonKarman(**kwargs, doDelta=True)
    np.testing.assert_almost_equal(vkd.kValue(1e10, 0).real, vkd.deltaAmplitude)

    # Either way, the fluxes should be the same.
    np.testing.assert_almost_equal(vk.flux, vkd.flux)
    assert vk != vkd
    assert vk.halfLightRadius != vkd.halfLightRadius


if __name__ == "__main__":
    test_vk()
    test_vk_delta()
