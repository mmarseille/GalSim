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
import os
import sys

import numpy as np
import galsim.deprecated.atmosphere

try:
    import galsim
except ImportError:
    path, filename = os.path.split(__file__)
    sys.path.append(os.path.abspath(os.path.join(path, "..")))
    import galsim

imgdir = os.path.join(".", "SBProfile_comparison_images") # Directory containing the reference
                                                          # images.

# AtmosphericPSF / Kolmogorov params and reference values
test_fwhm = 1.9
test_lor0 = 1.9
test_oversampling = 1.7

atmos_ref_fwhm_from_lor0 = test_lor0 * 0.976
atmos_ref_lor0_from_fwhm = test_fwhm / 0.976

# for flux normalization tests
test_flux = 1.9

# decimal point to go to for parameter value comparisons
param_decimal = 12

def funcname():
    import inspect
    return inspect.stack()[1][3]


def test_AtmosphericPSF_properties():
    """Test some basic properties of a known Atmospheric PSF.
    """
    import time
    t1 = time.time()
    apsf = galsim.deprecated.AtmosphericPSF(lam_over_r0=1.5)
    # Check that we are centered on (0, 0)
    cen = galsim._galsim.PositionD(0, 0)
    np.testing.assert_array_almost_equal(
            [apsf.centroid().x, apsf.centroid().y], [cen.x, cen.y], 10,
            err_msg="Atmospheric PSF not centered on (0, 0)")
    # Check Fourier properties
    np.testing.assert_almost_equal(apsf.maxK(), 5.8341564391716183, 9,
                                   err_msg="Atmospheric PSF .maxk() does not return known value.")
    np.testing.assert_almost_equal(apsf.stepK(), 1.0275679547331542, 9,
                                   err_msg="Atmospheric PSF .stepk() does not return known value.")
    np.testing.assert_almost_equal(apsf.kValue(cen), 1+0j, 4,
                                   err_msg="Atmospheric PSF k value at (0, 0) is not 1+0j.")
    t2 = time.time()
    print 'time for %s = %.2f'%(funcname(),t2-t1)

def test_AtmosphericPSF_flux():
    """Test that the flux of the atmospheric PSF is normalized to unity.
    """
    import time
    t1 = time.time()
    lors = np.linspace(0.5, 2., 5) # Different lambda_over_r0 values
    for lor in lors:
        apsf = galsim.deprecated.AtmosphericPSF(lam_over_r0=lor)
        print 'apsf.getFlux = ',apsf.getFlux()
        np.testing.assert_almost_equal(apsf.getFlux(), 1., 6,
                                       err_msg="Flux of atmospheric PSF (ImageViewD) is not 1.")
        # .draw() throws a warning if it doesn't get a float. This includes np.float64. Convert to
        # have the test pass.
        dx = float(lor / 10.)
        img = galsim.ImageF(256,256)
        img_array = apsf.draw(image=img, dx=dx).array
        np.testing.assert_almost_equal(img_array.sum(), 1., 3,
                                       err_msg="Flux of atmospheric PSF (image array) is not 1.")
    t2 = time.time()
    print 'time for %s = %.2f'%(funcname(),t2-t1)

def test_AtmosphericPSF_fwhm():
    """Test that the FWHM of the atmospheric PSF corresponds to the one expected from the
    lambda / r0 input."""
    import time
    t1 = time.time()
    lors = np.linspace(0.5, 2., 5) # Different lambda_over_r0 values
    for lor in lors:
        apsf = galsim.deprecated.AtmosphericPSF(lam_over_r0=lor)
        # .draw() throws a warning if it doesn't get a float. This includes np.float64. Convert to
        # have the test pass.
        dx_scale = 10
        dx = float(lor / dx_scale)
        # Need use_true_center=False, since we want the maximum to actually be drawn in one
        # of the pixels, rather than between the central 4 pixels.
        psf_array = apsf.draw(dx=dx, use_true_center=False).array
        nx, ny = psf_array.shape
        profile = psf_array[nx / 2, ny / 2:]
        # Now get the last array index where the profile value exceeds half the peak value as a
        # rough estimator of the HWHM.
        hwhm_index = np.where(profile > profile.max() / 2.)[0][-1]
        np.testing.assert_equal(hwhm_index, dx_scale / 2,
                                err_msg="Kolmogorov PSF does not have the expected FWHM.")
    t2 = time.time()
    print 'time for %s = %.2f'%(funcname(),t2-t1)

def test_atmos_flux_scaling():
    """Test flux scaling for AtmosphericPSF.
    """
    import time
    t1 = time.time()
    # init with lam_over_r0 and flux only (should be ok given last tests)
    obj = galsim.deprecated.AtmosphericPSF(lam_over_r0=test_lor0, flux=test_flux)
    obj *= 2.
    np.testing.assert_almost_equal(
        obj.getFlux(), test_flux * 2., decimal=param_decimal,
        err_msg="Flux param inconsistent after __imul__.")
    obj = galsim.deprecated.AtmosphericPSF(lam_over_r0=test_lor0, flux=test_flux)
    obj /= 2.
    np.testing.assert_almost_equal(
        obj.getFlux(), test_flux / 2., decimal=param_decimal,
        err_msg="Flux param inconsistent after __idiv__.")
    obj = galsim.deprecated.AtmosphericPSF(lam_over_r0=test_lor0, flux=test_flux)
    obj2 = obj * 2.
    # First test that original obj is unharmed...
    np.testing.assert_almost_equal(
        obj.getFlux(), test_flux, decimal=param_decimal,
        err_msg="Flux param inconsistent after __rmul__ (original).")
    # Then test new obj2 flux
    np.testing.assert_almost_equal(
        obj2.getFlux(), test_flux * 2., decimal=param_decimal,
        err_msg="Flux param inconsistent after __rmul__ (result).")
    obj = galsim.deprecated.AtmosphericPSF(lam_over_r0=test_lor0, flux=test_flux)
    obj2 = 2. * obj
    # First test that original obj is unharmed...
    np.testing.assert_almost_equal(
        obj.getFlux(), test_flux, decimal=param_decimal,
        err_msg="Flux param inconsistent after __mul__ (original).")
    # Then test new obj2 flux
    np.testing.assert_almost_equal(
        obj2.getFlux(), test_flux * 2., decimal=param_decimal,
        err_msg="Flux param inconsistent after __mul__ (result).")
    obj = galsim.deprecated.AtmosphericPSF(lam_over_r0=test_lor0, flux=test_flux)
    obj2 = obj / 2.
    # First test that original obj is unharmed...
    np.testing.assert_almost_equal(
        obj.getFlux(), test_flux, decimal=param_decimal,
        err_msg="Flux param inconsistent after __div__ (original).")
    # Then test new obj2 flux
    np.testing.assert_almost_equal(
        obj2.getFlux(), test_flux / 2., decimal=param_decimal,
        err_msg="Flux param inconsistent after __div__ (result).")
    t2 = time.time()
    print 'time for %s = %.2f'%(funcname(),t2-t1)


if __name__ == "__main__":
    test_AtmosphericPSF_flux()
    test_AtmosphericPSF_properties()
    test_AtmosphericPSF_fwhm()
    test_atmos_flux_scaling()
