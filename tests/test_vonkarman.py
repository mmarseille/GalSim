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


def test_vk_eq_kolm():
    lam = 500.0
    r0 = 0.2
    L0 = 1e3
    vk = galsim.VonKarman(lam, r0, L0)
    print(vk.lam)
    print(vk.r0)
    print(vk.L0)
    print(vk.kcrit)
    print(vk.stepK())
    print(vk.maxK())
    print(vk.gsparams)
    print(vk.structureFunction(0.0))
    print(vk.structureFunction(1e6))
    print(vk.kValue(0,1e12))

    vk.drawImage()

if __name__ == "__main__":
    test_vk_eq_kolm()
