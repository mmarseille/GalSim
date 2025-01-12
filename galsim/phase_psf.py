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
"""@file phase_psf.py
Utilities for creating PSFs from phase screens.

For PSFs drawn using real-space or Fourier methods, these utilities essentially evaluate the Fourier
optics diffraction equation:

PSF(x, y) = int( |FT(aperture(u, v) * exp(i * phase(u, v, x, y, t)))|^2, dt)

where x, y are focal plane coordinates and u, v are pupil plane coordinates.

For method='phot', one of two possible strategies are available.  The first strategy is to draw the
PSF using Fourier methods into an InterpolatedImage, and then shoot photons from that profile.  This
strategy has good accuracy, but can be computationally expensive, particularly for atmospheric PSFs
that need to be built up in small increments to simulate a finite exposure time. The second
strategy, which can be significantly faster, especially for atmospheric PSFs, is to use the
geometric optics approximation.  This approximation has good accuracy for atmospheric PSFs, so we
make it the default for PhaseScreenPSF.  The accuracy is somewhat less good for purely optical PSFs
though, so the default behavior for OpticalPSF is to use the first strategy.  The
`geometric_shooting` keyword can be used in both cases to override the default.


The main classes of note are:

Aperture
  Class representing the illuminated region of pupil.

AtmosphericScreen
  Class implementing phase(u, v, x, y, t) for von Karman type turbulence, with possibly evolving
  "non-frozen-flow" phases.

OpticalScreen
  Class implementing optical aberrations using Zernike polynomial expansions in the wavefront.

PhaseScreenList
  Python sequence type to hold multiple phase screens, for instance to simulate turbulence at
  different altitudes, or self-consistently model atmospheric and optical phase aberrations.  A key
  method is makePSF(), which will take the list of phase screens, add them together linearly
  (Fraunhofer approximation), and evaluate the above diffraction equation to yield a
  PhaseScreenPSF object.

PhaseScreenPSF
  A GSObject holding the evaluated PSF from a set of phase screens.

OpticalPSF
  A GSObject for optical PSFs with potentially complicated pupils and Zernike aberrations.

Atmosphere
  Convenience function to quickly assemble multiple AtmosphericScreens into a PhaseScreenList.
"""

from itertools import chain
from builtins import range, str
from heapq import heappush, heappop

import numpy as np
import galsim
from galsim import GSObject


class Aperture(object):
    """ Class representing a telescope aperture embedded in a larger pupil plane array -- for use
    with the PhaseScreenPSF class to create PSFs via Fourier or geometric optics.

    The pupil plane array is completely specified by its size, sampling interval, and pattern of
    illuminated pixels.  Pupil plane arrays can be specified either geometrically or using an image
    to indicate the illuminated pixels.  In both cases, various options exist to control the pupil
    plane size and sampling interval.

    Geometric pupil specification
    -----------------------------

    The first way to specify the details of the telescope aperture is through a series of keywords
    indicating the diameter, size of the central obscuration, and the nature of the struts
    holding up the secondary mirror (or prime focus cage, etc.).  The struts are assumed to be
    rectangular obscurations extending from the outer edge of the pupil to the outer edge of the
    obscuration disk (or to the pupil center if `obscuration = 0.`).  You can specify how many
    struts there are (evenly spaced in angle), how thick they are as a fraction of the pupil
    diameter, and what angle they start at relative to the positive y direction.

    The size (in meters) and sampling interval (in meters) of the pupil plane array representing the
    aperture can be set directly using the the `pupil_plane_size` and `pupil_plane_scale` keywords.
    However, in most situations, it's probably more convenient to let GalSim set these automatically
    based on the pupil geometry and the nature of the (potentially time-varying) phase aberrations
    from which a PSF is being derived.

    The pupil plane array physical size is by default set to twice the pupil diameter producing a
    Nyquist sampled PSF image.  While this would always be sufficient if using sinc interpolation
    over the PSF image for subsequent operations, GalSim by default uses the much faster (though
    approximate) quintic interpolant, which means that in some cases -- in particular, for
    significantly aberrated optical PSFs without atmospheric aberrations -- it may be useful to
    further increase the size of the pupil plane array, thereby increasing the sampling rate of the
    resulting PSF image.  This can be done by increasing the `oversampling` keyword.

    A caveat to the above occurs when using `geometric_shooting=True` to draw using photon-shooting.
    In this case, we only need an array just large enough to avoid clipping the pupil, which we can
    get by setting `oversampling=0.5`.

    The pupil plane array physical sampling interval (which is directly related to the resulting PSF
    image physical size) is set by default to the same interval as would be used to avoid
    significant aliasing (image folding) for an obscured Airy profile with matching diameter and
    obscuration and for the value of `folding_threshold` in the optionally specified gsparams
    argument.  If the phase aberrations are significant, however, the PSF image size computed this
    way may still not be sufficiently large to avoid aliasing.  To further increase the pupil plane
    sampling rate (and hence the PSF image size), you can increase the value of the `pad_factor`
    keyword.  An additional way to set the pupil sampling interval for a particular set of phase
    screens (i.e., for a particular PhaseScreenList) is to provide the screens in the `screen_list`
    argument.  Each screen in the list computes its own preferred sampling rate and the
    PhaseScreenList appropriately aggregates these. This last option also requires that a wavelength
    `lam` be specified, and is particularly helpful for creating PSFs derived from turbulent
    atmospheric screens.

    Finally, when specifying the pupil geometrically, Aperture may choose to make a small adjustment
    to `pupil_plane_scale` in order to produce an array with a good size for FFTs.  If your
    application depends on knowing the size and scale used with the Fourier optics framework, you
    can obtain these from the `aper.pupil_plane_size` and `aper.pupil_plane_scale` attributes.

    Pupil image specification
    --------------------------

    The second way to specify the pupil plane configuration is by passing in an image of it.  This
    can be useful, for example, if the struts are not evenly spaced or are not radially directed, as
    is assumed by the simple model for struts described above.  In this case, an exception is raised
    if keywords related to struts are also given.  On the other hand, the `obscuration` keyword is
    still used to ensure that the PSF images are not aliased, though it is ignored during the actual
    construction of the pupil plane illumination pattern.  Note that for complicated pupil
    configurations, it may be desireable to increase `pad_factor` for more fidelity at the expense
    of slower running time.  Finally, the `pupil_plane_im` that is passed in can be rotated during
    internal calculations by specifying a `pupil_angle` keyword.

    If you choose to pass in a pupil plane image, it must be a square array in which the image of
    the pupil is centered.  The areas that are illuminated should have some value >0, and the other
    areas should have a value of precisely zero.  Based on what the Aperture class determines is a
    good PSF sampling interval, the image of the pupil plane that is passed in might be zero-padded
    during internal calculations.  (The pupil plane array size and scale values can be accessed via
    the `aper.pupil_plane_size` and `aper.pupil_plane_scale` attributes.) The pixel scale of
    the pupil plane can be specified in one of three ways.  In descending order of priority, these
    are:
      1.  The `pupil_plane_scale` keyword argument (units are meters).
      2.  The `pupil_plane_im.scale` attribute (units are meters).
      3.  If (1) and (2) are both None, then the scale will be inferred by assuming that the
          illuminated pixel farthest from the image center is at a physical distance of self.diam/2.

    The `pupil_plane_size` and `lam` keywords are both ignored when constructing an Aperture from an
    image.

    @param diam                Aperture diameter in meters.
    @param lam                 Wavelength in nanometers.  [default: None]
    @param circular_pupil      Adopt a circular pupil?  [default: True]
    @param obscuration         Linear dimension of central obscuration as fraction of aperture
                               linear dimension. [0., 1.).  [default: 0.0]
    @param nstruts             Number of radial support struts to add to the central obscuration.
                               [default: 0]
    @param strut_thick         Thickness of support struts as a fraction of aperture diameter.
                               [default: 0.05]
    @param strut_angle         Angle made between the vertical and the strut starting closest to it,
                               defined to be positive in the counter-clockwise direction; must be an
                               Angle instance. [default: 0. * galsim.degrees]
    @param oversampling        Optional oversampling factor *in the image plane* for the PSF
                               eventually constructed using this Aperture.  Setting
                               `oversampling < 1` will produce aliasing in the PSF (not good).
                               [default: 1.0]
    @param pad_factor          Additional multiple by which to extend the PSF image to avoid
                               folding.  [default: 1.0]
    @param screen_list         An optional PhaseScreenList object.  If present, then get a good
                               pupil sampling interval using this object.  [default: None]
    @param pupil_plane_im      The GalSim.Image, NumPy array, or name of file containing the pupil
                               plane image, to be used instead of generating one based on the
                               obscuration and strut parameters.  [default: None]
    @param pupil_angle         If `pupil_plane_im` is not None, rotation angle for the pupil plane
                               (positive in the counter-clockwise direction).  Must be an Angle
                               instance. [default: 0. * galsim.degrees]
    @param pupil_plane_scale   Sampling interval in meters to use for the pupil plane array.  In
                               most cases, it's a good idea to leave this as None, in which case
                               GalSim will attempt to find a good value automatically.  The
                               exception is when specifying the pupil arrangement via an image, in
                               which case this keyword can be used to indicate the sampling of that
                               image.  See also `pad_factor` for adjusting the pupil sampling scale.
                               [default: None]
    @param pupil_plane_size    Size in meters to use for the pupil plane array.  In most cases, it's
                               a good idea to leave this as None, in which case GalSim will attempt
                               to find a good value automatically.  See also `oversampling` for
                               adjusting the pupil size.  [default: None]
    @param gsparams            An optional GSParams argument.  See the docstring for GSParams for
                               details. [default: None]
    """
    def __init__(self, diam, lam=None, circular_pupil=True, obscuration=0.0,
                 nstruts=0, strut_thick=0.05, strut_angle=0.0*galsim.degrees,
                 oversampling=1.0, pad_factor=1.0, screen_list=None,
                 pupil_plane_im=None, pupil_angle=0.0*galsim.degrees,
                 pupil_plane_scale=None, pupil_plane_size=None,
                 gsparams=None):

        self.diam = diam  # Always need to explicitly specify an aperture diameter.
        self._gsparams = gsparams

        if obscuration >= 1.:
            raise ValueError("Pupil fully obscured! obscuration = {:0} (>= 1)"
                             .format(obscuration))

        # You can either set geometric properties, or use a pupil image, but not both, so check for
        # that here.  One caveat is that we allow sanity checking the sampling of a pupil_image by
        # comparing it to the sampling GalSim would have used for an (obscured) Airy profile.  So
        # it's okay to specify an obscuration and a pupil_plane_im together, for example, but not
        # a pupil_plane_im and struts.
        is_default_geom = (circular_pupil and
                           nstruts == 0 and
                           strut_thick == 0.05 and
                           strut_angle == 0.0*galsim.degrees)
        if not is_default_geom and pupil_plane_im is not None:
            raise ValueError("Can't specify both geometric parameters and pupil_plane_im.")

        if screen_list is not None and lam is None:
            raise ValueError("Wavelength `lam` must be specified with `screen_list`.")

        # Although the user can set the pupil plane size and scale directly if desired, in most
        # cases it's nicer to have GalSim try to pick good values for these.

        # For the pupil plane size, we'll achieve Nyquist sampling in the focal plane if we sample
        # out to twice the diameter of the actual aperture in the pupil plane (completely
        # independent of wavelength, struts, obscurations, GSparams, and so on!).  This corresponds
        # to oversampling=1.0.  In fact, if we were willing to always use sinc interpolation, there
        # would never be any reason to go beyond this.  In practice, we usually use a faster, but
        # less accurate, quintic interpolant, which means we can benefit from improved sampling
        # (oversampling > 1.0) in some cases, especially when we're *not* modeling an atmosphere
        # which would otherwise tend to damp contributions at large k.
        good_pupil_size = 2 * diam * oversampling

        # For the pupil plane sampling interval, details like the obscuration and GSParams *are*
        # important as they affect the amount of aliasing encountered.  (An Airy profile has an
        # infinite extent in real space, so it *always* aliases at some level, more so with an
        # obscuration than without.  The GSParams settings indicate how much aliasing we're
        # willing to tolerate, so it's required here.)  To pick a good sampling interval, we start
        # with the interval that would be used for an obscured Airy GSObject profile.  If the
        # `screen_list` argument was supplied, then we also check its .stepk propertry, which
        # aggregates a good sampling interval from all of the wrapped PhaseScreens, and keep the
        # smaller stepk.
        if lam is None:
            # For Airy, pupil_plane_scale is independent of wavelength.  We could build an Airy with
            # lam_over_diam=1.0 and then alter the `good_pupil_scale = ...` line below
            # appropriately, but it's easier to just arbitrarily set `lam=500` if it wasn't set.
            lam = 500.0
        airy = galsim.Airy(diam=diam, lam=lam, obscuration=obscuration, gsparams=self._gsparams)
        stepk = airy.stepk
        if screen_list is not None:
            screen_list = galsim.PhaseScreenList(screen_list)
            stepk = min(stepk,
                        screen_list._stepK(lam=lam, diam=diam, obscuration=obscuration,
                                          gsparams=self._gsparams))
        good_pupil_scale = (stepk * lam * 1.e-9 * (galsim.radians / galsim.arcsec)
                            / (2 * np.pi * pad_factor))

        # Now that we have good candidate sizes and scales, we load or generate the pupil plane
        # array.
        if pupil_plane_im is not None:  # Use image of pupil plane
            self._load_pupil_plane(pupil_plane_im, pupil_angle, pupil_plane_scale,
                                   good_pupil_scale, good_pupil_size)
        else:  # Use geometric parameters.
            if pupil_plane_scale is not None:
                # Check input scale and warn if looks suspicious.
                if pupil_plane_scale > good_pupil_scale:  # pragma: no cover
                    import warnings
                    ratio = good_pupil_scale / pupil_plane_scale
                    warnings.warn("Input pupil_plane_scale may be too large for good sampling.\n"
                                  "Consider decreasing pupil_plane_scale by a factor %f, and/or "
                                  "check PhaseScreenPSF outputs for signs of folding in real "
                                  "space."%(1./ratio))
            else:
                pupil_plane_scale = good_pupil_scale
            if pupil_plane_size is not None:
                # Check input size and warn if looks suspicious
                if pupil_plane_size < good_pupil_size:  # pragma: no cover
                    import warnings
                    ratio = good_pupil_size / pupil_plane_size
                    warnings.warn("Input pupil_plane_size may be too small for good focal-plane"
                                  "sampling.\n"
                                  "Consider increasing pupil_plane_size by a factor %f, and/or "
                                  "check PhaseScreenPSF outputs for signs of undersampling."%ratio)
            else:
                pupil_plane_size = good_pupil_size
            self._generate_pupil_plane(circular_pupil, obscuration,
                                       nstruts, strut_thick, strut_angle,
                                       pupil_plane_scale, pupil_plane_size)

        # Check FFT size
        if self._gsparams is not None:
            maximum_fft_size = self._gsparams.maximum_fft_size
        else:
            maximum_fft_size = galsim.GSParams().maximum_fft_size
        if self.npix > maximum_fft_size:
            raise RuntimeError("Created pupil plane array that is too large, {0} "
                               "If you can handle the large FFT, you may update "
                               "gsparams.maximum_fft_size".format(self.npix))

    def _generate_pupil_plane(self, circular_pupil, obscuration, nstruts, strut_thick, strut_angle,
                              pupil_plane_scale, pupil_plane_size):
        """ Create an array of illuminated pixels parameterically.
        """
        ratio = pupil_plane_size/pupil_plane_scale
        # Fudge a little to prevent good_fft_size() from turning 512.0001 into 768.
        ratio *= (1.0 - 1.0/2**14)
        self.npix = galsim.Image.good_fft_size(int(np.ceil(ratio)))
        self.pupil_plane_size = pupil_plane_size
        # Shrink scale such that size = scale * npix exactly.
        self.pupil_plane_scale = pupil_plane_size / self.npix
        # Save params for str/repr
        self._circular_pupil = circular_pupil
        self._obscuration = obscuration
        self._nstruts = nstruts
        self._strut_thick = strut_thick
        self._strut_angle = strut_angle

        radius = 0.5*self.diam
        if circular_pupil:
            self._illuminated = (self.rsqr < radius**2)
            if obscuration > 0.:
                self._illuminated *= self.rsqr >= (radius*obscuration)**2
        else:
            self._illuminated = (np.abs(self.u) < radius) & (np.abs(self.v) < radius)
            if obscuration > 0.:
                self._illuminated *= ((np.abs(self.u) >= radius*obscuration) *
                                      (np.abs(self.v) >= radius*obscuration))

        if nstruts > 0:
            if not isinstance(strut_angle, galsim.Angle):
                raise TypeError("Input kwarg strut_angle must be a galsim.Angle instance.")
            # Add the initial rotation if requested, converting to radians.
            rot_u, rot_v = self.u, self.v
            if strut_angle.rad != 0.:
                rot_u, rot_v = galsim.utilities.rotate_xy(rot_u, rot_v, -strut_angle)
            rotang = 360. * galsim.degrees / nstruts
            # Then loop through struts setting to zero the regions which lie under the strut
            for istrut in range(nstruts):
                rot_u, rot_v = galsim.utilities.rotate_xy(rot_u, rot_v, -rotang)
                self._illuminated *= ((np.abs(rot_u) >= radius * strut_thick) + (rot_v < 0.0))

    def _load_pupil_plane(self, pupil_plane_im, pupil_angle, pupil_plane_scale, good_pupil_scale,
                          good_pupil_size):
        """ Create an array of illuminated pixels with appropriate size and scale from an input
        image of the pupil.  The basic strategy is:

        1.  Read in array.
        2.  Determine the scale.
        3.  Pad the input array with zeros to meet the requested pupil size.
        4.  Check that the pupil plane sampling interval is at least as small as requested.
        5.  Optionally rotate pupil plane.
        """
        # Handle multiple types of input: NumPy array, galsim.Image, or string for filename with
        # image.
        if isinstance(pupil_plane_im, np.ndarray):
            # Make it into an image.
            pupil_plane_im = galsim.Image(pupil_plane_im)
        elif isinstance(pupil_plane_im, galsim.Image):
            # Make sure not to overwrite input image.
            pupil_plane_im = pupil_plane_im.copy()
        else:
            # Read in image of pupil plane from file.
            pupil_plane_im = galsim.fits.read(pupil_plane_im)
        # scale = pupil_plane_im.scale # Interpret as either the pixel scale in meters, or None.
        pp_arr = pupil_plane_im.array
        self.npix = pp_arr.shape[0]

        # Sanity checks
        if pupil_plane_im.array.shape[0] != pupil_plane_im.array.shape[1]:
            raise ValueError("We require square input pupil plane arrays!")
        if pupil_plane_im.array.shape[0] % 2 == 1:
            raise ValueError("Even-sized input arrays are required for the pupil plane!")

        # Set the scale, priority is:
        # 1.  pupil_plane_scale kwarg
        # 2.  image.scale if not None
        # 3.  Use diameter and farthest illuminated pixel.
        if pupil_plane_scale is not None:
            self.pupil_plane_scale = pupil_plane_scale
        elif pupil_plane_im.scale is not None:
            self.pupil_plane_scale = pupil_plane_im.scale
        else:
            # If self.pupil_plane_scale is not set yet, then figure it out from the distance
            # of the farthest illuminated pixel from the image center and the aperture diameter.
            # below is essentially np.linspace(-0.5, 0.5, self.npix)
            u = np.fft.fftshift(np.fft.fftfreq(self.npix))
            u, v = np.meshgrid(u, u)
            r = np.hypot(u, v)
            rmax_illum = np.max(r*(pupil_plane_im.array > 0))
            self.pupil_plane_scale = self.diam / (2.0 * rmax_illum * self.npix)
        self.pupil_plane_size = self.pupil_plane_scale * self.npix

        # Check the pupil plane size here and bump it up if necessary.
        if self.pupil_plane_size < good_pupil_size:
            new_npix = galsim.Image.good_fft_size(int(np.ceil(
                    good_pupil_size/self.pupil_plane_scale)))
            pad_width = (new_npix-self.npix)//2
            pp_arr = np.pad(pp_arr, [(pad_width, pad_width)]*2, mode='constant')
            self.npix = new_npix
            self.pupil_plane_size = self.pupil_plane_scale * self.npix

        # Check sampling interval and warn if it's not good enough.
        if self.pupil_plane_scale > good_pupil_scale:  # pragma: no cover
            import warnings
            ratio = self.pupil_plane_scale / good_pupil_scale
            warnings.warn("Input pupil plane image may not be sampled well enough!\n"
                          "Consider increasing sampling by a factor %f, and/or check "
                          "PhaseScreenPSF outputs for signs of folding in real space."%ratio)

        if pupil_angle.rad == 0.:
            self._illuminated = pp_arr.astype(bool)
        else:
            # Rotate the pupil plane image as required based on the `pupil_angle`, being careful to
            # ensure that the image is one of the allowed types.  We ignore the scale.
            b = galsim._BoundsI(1,self.npix,1,self.npix)
            im = galsim._Image(pp_arr, b, galsim.PixelScale(1.))
            int_im = galsim.InterpolatedImage(im, x_interpolant='linear',
                                              calculate_stepk=False, calculate_maxk=False)
            int_im = int_im.rotate(pupil_angle)
            new_im = galsim.Image(pp_arr.shape[1], pp_arr.shape[0])
            new_im = int_im.drawImage(image=new_im, scale=1., method='no_pixel')
            pp_arr = new_im.array
            # Restore hard edges that might have been lost during the interpolation.  To do this, we
            # check the maximum value of the entries.  Values after interpolation that are >half
            # that maximum value are kept as nonzero (True), but those that are <half the maximum
            # value are set to zero (False).
            max_pp_val = np.max(pp_arr)
            pp_arr[pp_arr < 0.5*max_pp_val] = 0.
            self._illuminated = pp_arr.astype(bool)

    # Used in Aperture.__str__ and OpticalPSF.__str__
    def _geometry_str(self):
        s = ""
        if not self._circular_pupil:
            s += ", circular_pupil=False"
        if self._obscuration != 0.0:
            s += ", obscuration=%s"%self._obscuration
        if self._nstruts != 0:
            s += ", nstruts=%s"%self._nstruts
            if self._strut_thick != 0.05:
                s += ", strut_thick=%s"%self._strut_thick
            if self._strut_angle != 0*galsim.degrees:
                s += ", strut_angle=%s"%self._strut_angle
        return s

    def __str__(self):
        s = "galsim.Aperture(diam=%r"%self.diam
        if hasattr(self, '_circular_pupil'):  # Pupil was created geometrically, so use that here.
            s += self._geometry_str()
        s += ")"
        return s

    # Used in Aperture.__repr__ and OpticalPSF.__repr__
    def _geometry_repr(self):
        s = ""
        if not self._circular_pupil:
            s += ", circular_pupil=False"
        if self._obscuration != 0.0:
            s += ", obscuration=%r"%self._obscuration
        if self._nstruts != 0:
            s += ", nstruts=%r"%self._nstruts
            if self._strut_thick != 0.05:
                s += ", strut_thick=%r"%self._strut_thick
            if self._strut_angle != 0*galsim.degrees:
                s += ", strut_angle=%r"%self._strut_angle
        return s

    def __repr__(self):
        s = "galsim.Aperture(diam=%r"%self.diam
        if hasattr(self, '_circular_pupil'):  # Pupil was created geometrically, so use that here.
            s += self._geometry_repr()
            s += ", pupil_plane_scale=%r"%self.pupil_plane_scale
            s += ", pupil_plane_size=%r"%self.pupil_plane_size
        else:  # Pupil was created from image, so use that instead.
            # It's slightly less annoying to see an enormous stream of zeros fly by than an enormous
            # stream of Falses, so convert to int16.
            tmp = self.illuminated.astype(np.int16).tolist()
            s += ", pupil_plane_im=array(%r"%tmp+", dtype='int16')"
            s += ", pupil_plane_scale=%r"%self.pupil_plane_scale
        if hasattr(self, '_gsparams') and self._gsparams is not None:
            s += ", gsparams=%r"%self._gsparams
        s += ")"
        return s

    def __eq__(self, other):
        return (isinstance(other, galsim.Aperture) and
                self.diam == other.diam and
                self.pupil_plane_scale == other.pupil_plane_scale and
                np.array_equal(self.illuminated, other.illuminated))

    def __hash__(self):
        # Cache since self.illuminated may be large.
        if not hasattr(self, '_hash'):
            self._hash = hash(("galsim.Aperture", self.diam, self.pupil_plane_scale))
            self._hash ^= hash(tuple(self.illuminated.ravel()))
        return self._hash

    # Properties show up nicely in the interactive terminal for
    #     >>>help(Aperture)
    # So we make a thin wrapper here.
    @property
    def illuminated(self):
        """  A boolean array indicating which positions in the pupil plane are exposed to the sky.
        """
        return self._illuminated

    @property
    def rho(self):
        """ Unit-disk normalized pupil plane coordinate as a complex number:
        (x, y) => x + 1j * y.
        """
        if not hasattr(self, '_rho') or self._rho is None:
            u = np.fft.fftshift(np.fft.fftfreq(self.npix, self.diam/self.pupil_plane_size/2.0))
            u, v = np.meshgrid(u, u)
            self._rho = u + 1j * v
        return self._rho

    @property
    def u(self):
        """Pupil horizontal coordinate array in meters."""
        if not hasattr(self, '_u'):
            u = np.fft.fftshift(np.fft.fftfreq(self.npix, 1./self.pupil_plane_size))
            self._u, self._v = np.meshgrid(u, u)
        return self._u

    @property
    def v(self):
        """Pupil vertical coordinate array in meters."""
        if not hasattr(self, '_v'):
            u = np.fft.fftshift(np.fft.fftfreq(self.npix, 1./self.pupil_plane_size))
            self._u, self._v = np.meshgrid(u, u)
        return self._v

    @property
    def rsqr(self):
        """Pupil radius squared array in meters squared."""
        if not hasattr(self, '_rsqr'):
            self._rsqr = self.u**2 + self.v**2
        return self._rsqr

    def __getstate__(self):
        # Let unpickled object reconstruct cached values on-the-fly instead of including them in the
        # pickle.
        d = self.__dict__
        for k in ['_rho', '_u', '_v', '_rsqr']:
            d.pop(k, None)
        return d

    # Some quick notes for Josh:
    # - Relation between real-space grid with size theta and pitch dtheta (dimensions of angle)
    #   and corresponding (fast) Fourier grid with size 2*maxk and pitch stepk (dimensions of
    #   inverse angle):
    #     stepk = 2*pi/theta
    #     maxk = pi/dtheta
    # - Relation between aperture of size L and pitch dL (dimensions of length, not angle!) and
    #   (fast) Fourier grid:
    #     dL = stepk * lambda / (2 * pi)
    #     L = maxk * lambda / pi
    # - Implies relation between aperture grid and real-space grid:
    #     dL = lambda/theta
    #     L = lambda/dtheta
    def _stepK(self, lam, scale_unit=galsim.arcsec):
        """Return the Fourier grid spacing for this aperture at given wavelength.

        @param lam         Wavelength in nanometers.
        @param scale_unit  Inverse units in which to return result [default: galsim.arcsec]
        @returns           Fourier grid spacing.
        """
        return 2*np.pi*self.pupil_plane_scale/(lam*1e-9) * scale_unit/galsim.radians

    def _maxK(self, lam, scale_unit=galsim.arcsec):
        """Return the Fourier grid half-size for this aperture at given wavelength.

        @param lam         Wavelength in nanometers.
        @param scale_unit  Inverse units in which to return result [default: galsim.arcsec]
        @returns           Fourier grid half-size.
        """
        return np.pi*self.pupil_plane_size/(lam*1e-9) * scale_unit/galsim.radians

    def _sky_scale(self, lam, scale_unit=galsim.arcsec):
        """Return the image scale for this aperture at given wavelength.
        @param lam         Wavelength in nanometers.
        @param scale_unit  Units in which to return result [default: galsim.arcsec]
        @returns           Image scale.
        """
        return (lam*1e-9) / self.pupil_plane_size * galsim.radians/scale_unit

    def _sky_size(self, lam, scale_unit=galsim.arcsec):
        """Return the image size for this aperture at given wavelength.
        @param lam         Wavelength in nanometers.
        @param scale_unit  Units in which to return result [default: galsim.arcsec]
        @returns           Image size.
        """
        return (lam*1e-9) / self.pupil_plane_scale * galsim.radians/scale_unit


class PhaseScreenList(object):
    """ List of phase screens that can be turned into a PSF.  Screens can be either atmospheric
    layers or optical phase screens.  Generally, one would assemble a PhaseScreenList object using
    the function `Atmosphere`.  Layers can be added, removed, appended, etc. just like items can be
    manipulated in a python list.  For example:

        # Create an atmosphere with three layers.
        >>> screens = galsim.PhaseScreenList([galsim.AtmosphericScreen(...),
                                              galsim.AtmosphericScreen(...),
                                              galsim.AtmosphericScreen(...)])
        # Add another layer
        >>> screens.append(galsim.AtmosphericScreen(...))
        # Remove the second layer
        >>> del screens[1]
        # Switch the first and second layer.  Silly, but works...
        >>> screens[0], screens[1] = screens[1], screens[0]

    Methods
    -------
    makePSF()          Obtain a PSF from this set of phase screens.  See PhaseScreenPSF docstring
                       for more details.
    wavefront()        Compute the cumulative wavefront due to all screens.
    wavefront_gradient()   Compute the cumulative wavefront gradient due to all screens.

    @param layers  Sequence of phase screens.
    """
    def __init__(self, *layers):
        if len(layers) == 1:
            # First check if layers[0] is a PhaseScreenList, so we avoid nesting.
            if isinstance(layers[0], galsim.PhaseScreenList):
                layers = layers[0]._layers
            else:
                # Next, see if layers[0] is iterable.  E.g., to catch generator expressions.
                try:
                    layers = list(layers[0])
                except TypeError:
                    # If that fails, check if layers[0] is a bare PhaseScreen.  Should probably
                    # make an ABC for this (use __subclasshook__?), but for now, just check
                    # AtmosphericScreen and OpticalScreen.
                    if isinstance(layers[0], (galsim.AtmosphericScreen, galsim.OpticalScreen)):
                        layers = [layers[0]]
        # else, layers is either empty or a tuple of PhaseScreens and so responds appropriately
        # to list() below.
        self._layers = list(layers)
        self._update_attrs()
        self._pending = []  # Pending PSFs to calculate upon first drawImage.
        self._update_time_heap = []  # Heap to store each PSF's next time-of-update.

    def __len__(self):
        return len(self._layers)

    def __getitem__(self, index):
        import numbers
        cls = type(self)
        if isinstance(index, slice):
            return cls(self._layers[index])
        elif isinstance(index, numbers.Integral):
            return self._layers[index]
        else:  # pragma: no cover
            msg = "{cls.__name__} indices must be integers"
            raise TypeError(msg.format(cls=cls))

    def __setitem__(self, index, layer):
        self._layers[index] = layer
        self._update_attrs()

    def __delitem__(self, index):
        del self._layers[index]
        self._update_attrs()

    def append(self, layer):
        self._layers.append(layer)
        self._update_attrs()

    def extend(self, layers):
        self._layers.extend(layers)
        self._update_attrs()

    def __str__(self):
        return "galsim.PhaseScreenList([%s])" % ",".join(str(l) for l in self._layers)

    def __repr__(self):
        return "galsim.PhaseScreenList(%r)" % self._layers

    def __eq__(self, other):
        return isinstance(other,PhaseScreenList) and self._layers == other._layers

    def __ne__(self, other): return not self == other

    __hash__ = None  # Mutable means not hashable.

    def _update_attrs(self):
        # If any of the wrapped PhaseScreens have an rng, then eval(repr(screen_list)) will run, but
        # fail to round-trip to the original object.  So we search for that here and set/delete a
        # dummy rng sentinel attribute so do_pickle() will know to skip the obj == eval(repr(obj))
        # test.
        self.__dict__.pop('rng', None)
        if any(hasattr(l, 'rng') for l in self):
            self.rng = None
        self.dynamic = any(l.dynamic for l in self)
        self.reversible = all(l.reversible for l in self)

    def _seek(self, t):
        """Set all layers' internal clocks to time t."""
        for layer in self:
            try:
                layer._seek(t)
            except AttributeError:
                # Time indep phase screen
                pass
        self._update_attrs()

    def _reset(self):
        """Reset phase screens back to time=0."""
        for layer in self:
            try:
                layer._reset()
            except AttributeError:
                # Time indep phase screen
                pass
        self._update_attrs()

    def _delayCalculation(self, psf):
        """Add psf to delayed calculation list."""
        self._pending.append(psf)
        heappush(self._update_time_heap, (psf.t0, len(self._pending)-1))

    def _prepareDraw(self):
        """Calculate previously delayed PSFs."""
        if not self._pending:
            return
            # See if we have any dynamic screens.  If not, then we can immediately compute each PSF
            # in a simple loop.
        if not self.dynamic:
            for psf in self._pending:
                psf._step()
                psf._finalize()
            self._pending = []
            self._update_time_heap = []
            return

        # If we do have time-evolving screens, then iteratively increment the time while being
        # careful to always stop at multiples of each PSF's time_step attribute to update that PSF.
        # Use a heap to track the next time to stop at.
        while(self._update_time_heap):
            # Get and seek to next time that has a PSF update.
            t, i = heappop(self._update_time_heap)
            self._seek(t)
            # Update that PSF
            psf = self._pending[i]
            psf._step()
            # If that PSF's next possible update time doesn't extend past its exptime, then
            # push it back on the heap.
            t += psf.time_step
            if t < psf.t0 + psf.exptime:
                heappush(self._update_time_heap, (t, i))
            else:
                psf._finalize()
        self._pending = []

    def wavefront(self, u, v, t, theta=(0.0*galsim.arcmin, 0.0*galsim.arcmin)):
        """ Compute cumulative wavefront due to all phase screens in PhaseScreenList.

        Wavefront here indicates the distance by which the physical wavefront lags or leads the
        ideal plane wave (pre-optics) or spherical wave (post-optics).

        @param u        Horizontal pupil coordinate (in meters) at which to evaluate wavefront.  Can
                        be a scalar or an iterable.  The shapes of u and v must match.
        @param v        Vertical pupil coordinate (in meters) at which to evaluate wavefront.  Can
                        be a scalar or an iterable.  The shapes of u and v must match.
        @param t        Times (in seconds) at which to evaluate wavefront.  Can be a scalar or an
                        iterable.  If scalar, then the size will be broadcast up to match that of
                        u and v.  If iterable, then the shape must match the shapes of u and v.
        @param theta    Field angle at which to evaluate wavefront, as a 2-tuple of `galsim.Angle`s.
                        [default: (0.0*galsim.arcmin, 0.0*galsim.arcmin)]  Only a single theta is
                        permitted.
        @returns        Array of wavefront lag or lead in nanometers.
        """
        if len(self._layers) > 1:
            return np.sum([layer.wavefront(u, v, t, theta) for layer in self], axis=0)
        else:
            return self._layers[0].wavefront(u, v, t, theta)

    def wavefront_gradient(self, u, v, t, theta=(0.0*galsim.arcmin, 0.0*galsim.arcmin)):
        """ Compute cumulative wavefront gradient due to all phase screens in PhaseScreenList.

        @param u        Horizontal pupil coordinate (in meters) at which to evaluate wavefront.  Can
                        be a scalar or an iterable.  The shapes of u and v must match.
        @param v        Vertical pupil coordinate (in meters) at which to evaluate wavefront.  Can
                        be a scalar or an iterable.  The shapes of u and v must match.
        @param t        Times (in seconds) at which to evaluate wavefront.  Can be a scalar or an
                        iterable.  If scalar, then the size will be broadcast up to match that of
                        u and v.  If iterable, then the shape must match the shapes of u and v.
        @param theta    Field angle at which to evaluate wavefront, as a 2-tuple of `galsim.Angle`s.
                        [default: (0.0*galsim.arcmin, 0.0*galsim.arcmin)]  Only a single theta is
                        permitted.
        @returns        Arrays dWdu and dWdv of wavefront lag or lead gradient in nm/m.
        """
        if len(self._layers) > 1:
            return np.sum([layer.wavefront_gradient(u, v, t, theta) for layer in self], axis=0)
        else:
            return self._layers[0].wavefront_gradient(u, v, t, theta)

    def _wavefront(self, u, v, t, theta):
        if len(self._layers) > 1:
            return np.sum([layer._wavefront(u, v, t, theta) for layer in self], axis=0)
        else:
            return self._layers[0]._wavefront(u, v, t, theta)

    def _wavefront_gradient(self, u, v, t, theta):
        if len(self._layers) > 1:
            return np.sum([layer._wavefront_gradient(u, v, t, theta) for layer in self], axis=0)
        else:
            return self._layers[0]._wavefront_gradient(u, v, t, theta)

    def makePSF(self, lam, **kwargs):
        """Create a PSF from the current PhaseScreenList.

        @param lam                 Wavelength in nanometers at which to compute PSF.
        @param t0                  Time at which to start exposure in seconds.  [default: 0.0]
        @param exptime             Time in seconds over which to accumulate evolving instantaneous
                                   PSF.  [default: 0.0]
        @param time_step           Time interval in seconds with which to sample phase screens when
                                   drawing using real-space or Fourier methods, or when using
                                   photon-shooting without the geometric optics approximation.  Note
                                   that the default value of 0.025 is fairly arbitrary.  For careful
                                   studies, we recommend checking that results are stable when
                                   decreasing time_step.  Also note that when drawing using
                                   photon-shooting with the geometric optics approximation this
                                   keyword is ignored, as the phase screen can be sampled
                                   continuously in this case instead of at discrete intervals.
                                   [default: 0.025]
        @param flux                Flux of output PSF.  [default: 1.0]
        @param theta               Field angle of PSF as a 2-tuple of Angles.
                                   [default: (0.0*galsim.arcmin, 0.0*galsim.arcmin)]
        @param interpolant         Either an Interpolant instance or a string indicating which
                                   interpolant should be used.  Options are 'nearest', 'sinc',
                                   'linear', 'cubic', 'quintic', or 'lanczosN' where N should be the
                                   integer order to use. [default: galsim.Quintic()]
        @param scale_unit          Units to use for the sky coordinates of the output profile.
                                   [default: galsim.arcsec]
        @param ii_pad_factor       Zero-padding factor by which to extend the image of the PSF when
                                   creating the `InterpolatedImage`.  See the `InterpolatedImage`
                                   docstring for more details.  [default: 4.]
        @param suppress_warning    If `pad_factor` is too small, the code will emit a warning
                                   telling you its best guess about how high you might want to raise
                                   it.  However, you can suppress this warning by using
                                   `suppress_warning=True`.  [default: False]
        @param geometric_shooting  If True, then when drawing using photon shooting, use geometric
                                   optics approximation where the photon angles are derived from the
                                   phase screen gradient.  If False, then first draw using Fourier
                                   optics and then shoot from the derived InterpolatedImage.
                                   [default: True]
        @param aper                Aperture to use to compute PSF(s).  [default: None]
        @param gsparams            An optional GSParams argument.  See the docstring for GSParams
                                   for details.  [default: None]

        The following are optional keywords to use to setup the aperture if `aper` is not provided.

        @param diam                Aperture diameter in meters.
        @param circular_pupil      Adopt a circular pupil?  [default: True]
        @param obscuration         Linear dimension of central obscuration as fraction of aperture
                                   linear dimension. [0., 1.).  [default: 0.0]
        @param nstruts             Number of radial support struts to add to the central
                                   obscuration. [default: 0]
        @param strut_thick         Thickness of support struts as a fraction of aperture diameter.
                                   [default: 0.05]
        @param strut_angle         Angle made between the vertical and the strut starting closest to
                                   it, defined to be positive in the counter-clockwise direction;
                                   must be an Angle instance. [default: 0. * galsim.degrees]
        @param oversampling        Optional oversampling factor *in the image plane* for the PSF
                                   eventually constructed using this Aperture.  Setting
                                   `oversampling < 1` will produce aliasing in the PSF (not good).
                                   [default: 1.0]
        @param pad_factor          Additional multiple by which to extend the PSF image to avoid
                                   folding.  [default: 1.0]
        @param pupil_plane_im      The GalSim.Image, NumPy array, or name of file containing the
                                   pupil plane image, to be used instead of generating one based on
                                   the obscuration and strut parameters.  [default: None]
        @param pupil_angle         If `pupil_plane_im` is not None, rotation angle for the pupil
                                   plane (positive in the counter-clockwise direction).  Must be an
                                   Angle instance. [default: 0. * galsim.degrees]
        @param pupil_plane_scale   Sampling interval in meters to use for the pupil plane array.  In
                                   most cases, it's a good idea to leave this as None, in which case
                                   GalSim will attempt to find a good value automatically.  The
                                   exception is when specifying the pupil arrangement via an image,
                                   in which case this keyword can be used to indicate the sampling
                                   of that image.  See also `pad_factor` for adjusting the pupil
                                   sampling scale. [default: None]
        @param pupil_plane_size    Size in meters to use for the pupil plane array.  In most cases,
                                   it's a good idea to leave this as None, in which case GalSim will
                                   attempt to find a good value automatically.  See also
                                   `oversampling` for adjusting the pupil size.  [default: None]
        """
        # Determine if theta is a single 2-tuple of Angles (okay) or an iterable of 2-tuples of
        # Angles (deprecated).
        theta = kwargs.pop('theta', (0.0*galsim.arcmin, 0.0*galsim.arcmin))

        # 2-tuples are iterable, so to check whether theta is indicating a single pointing, or a
        # generator of pointings we need to look at the first item.  If the first item is
        # iterable itself, then assume theta is an iterable of 2-tuple field angles.  We then
        # replace the consumed tuple at the beginning of the generator and go on.  If the first
        # item is scalar, then assume that it's the x-component of a single field angle.
        theta = iter(theta)
        th0 = next(theta)
        if not hasattr(th0, '__iter__'):
            theta = [th0, next(theta)]
            return PhaseScreenPSF(self, lam, theta=theta, **kwargs)
        else:
            from .deprecated import depr
            depr('list of `theta`s', 1.5, '[psl.makePSF(..., theta=th) for th in theta]')
            theta = chain([th0], theta)
            return [PhaseScreenPSF(self, lam, theta=th, **kwargs) for th in theta]

    @property
    def r0_500_effective(self):
        """Effective r0_500 for set of screens in list that define an r0_500 attribute."""
        return np.sum([l.r0_500**(-5./3) for l in self if hasattr(l, 'r0_500')])**(-3./5)

    def _stepK(self, **kwargs):
        """Return an appropriate stepk for this list of phase screens.

        The required set of parameters depends on the types of the individual PhaseScreens in the
        PhaseScreenList.  See the documentation for the individual PhaseScreen.pupil_plane_scale
        methods for more details.

        @returns  stepk.
        """
        # Generically, GalSim propagates stepk for convolutions using
        #   stepk = sum(s**-2 for s in stepks)**(-0.5)
        # We're not actually doing convolution between screens here, though.  In fact, the right
        # relation for Kolmogorov screens uses exponents -5./3 and -3./5:
        #   stepk = sum(s**(-5./3) for s in stepks)**(-3./5)
        # Since most of the layers in a PhaseScreenList are likely to be (nearly) Kolmogorov
        # screens, we'll use that relation.
        return np.sum([layer._stepK(**kwargs)**(-5./3) for layer in self])**(-3./5)


class PhaseScreenPSF(GSObject):
    """A PSF surface brightness profile constructed by integrating over time the instantaneous PSF
    derived from a set of phase screens and an aperture.

    There are two equivalent ways to construct a PhaseScreenPSF given a PhaseScreenList:
        >>> psf = screen_list.makePSF(...)
        >>> psf = PhaseScreenPSF(screen_list, ...)

    Computing a PSF from a phase screen also requires an Aperture be specified.  This can be done
    either directly via the `aper` keyword, or by setting a number of keywords that will be passed
    to the `Aperture` constructor.  The `aper` keyword always takes precedence.

    @param screen_list         PhaseScreenList object from which to create PSF.
    @param lam                 Wavelength in nanometers at which to compute PSF.
    @param t0                  Time at which to start exposure in seconds.  [default: 0.0]
    @param exptime             Time in seconds over which to accumulate evolving instantaneous PSF.
                               [default: 0.0]
    @param time_step           Time interval in seconds with which to sample phase screens when
                               drawing using real-space or Fourier methods, or when using
                               photon-shooting without the geometric optics approximation.  Note
                               that the default value of 0.025 is fairly arbitrary.  For careful
                               studies, we recommend checking that results are stable when
                               decreasing time_step.  Also note that when drawing using
                               photon-shooting with the geometric optics approximation this
                               keyword is ignored, as the phase screen can be sampled
                               continuously in this case instead of at discrete intervals.
                               [default: 0.025]
    @param flux                Flux of output PSF [default: 1.0]
    @param theta               Field angle of PSF as a 2-tuple of Angles.
                               [default: (0.0*galsim.arcmin, 0.0*galsim.arcmin)]
    @param interpolant         Either an Interpolant instance or a string indicating which
                               interpolant should be used.  Options are 'nearest', 'sinc', 'linear',
                               'cubic', 'quintic', or 'lanczosN' where N should be the integer order
                               to use.  [default: galsim.Quintic()]
    @param scale_unit          Units to use for the sky coordinates of the output profile.
                               [default: galsim.arcsec]
    @param ii_pad_factor       Zero-padding factor by which to extend the image of the PSF when
                               creating the `InterpolatedImage`.  See the `InterpolatedImage`
                               docstring for more details.  [default: 4.]
    @param suppress_warning    If `pad_factor` is too small, the code will emit a warning telling
                               you its best guess about how high you might want to raise it.
                               However, you can suppress this warning by using
                               `suppress_warning=True`.  [default: False]
    @param geometric_shooting  If True, then when drawing using photon shooting, use geometric
                               optics approximation where the photon angles are derived from the
                               phase screen gradient.  If False, then first draw using Fourier
                               optics and then shoot from the derived InterpolatedImage.
                               [default: True]
    @param aper                Aperture to use to compute PSF(s).  [default: None]
    @param gsparams            An optional GSParams argument.  See the docstring for GSParams for
                               details. [default: None]

    The following are optional keywords to use to setup the aperture if `aper` is not provided:

    @param diam                Aperture diameter in meters.
    @param circular_pupil      Adopt a circular pupil?  [default: True]
    @param obscuration         Linear dimension of central obscuration as fraction of aperture
                               linear dimension. [0., 1.).  [default: 0.0]
    @param nstruts             Number of radial support struts to add to the central obscuration.
                               [default: 0]
    @param strut_thick         Thickness of support struts as a fraction of aperture diameter.
                               [default: 0.05]
    @param strut_angle         Angle made between the vertical and the strut starting closest to it,
                               defined to be positive in the counter-clockwise direction; must be an
                               Angle instance. [default: 0. * galsim.degrees]
    @param oversampling        Optional oversampling factor *in the image plane* for the PSF
                               eventually constructed using this Aperture.  Setting
                               `oversampling < 1` will produce aliasing in the PSF (not good).
                               [default: 1.0]
    @param pad_factor          Additional multiple by which to extend the PSF image to avoid
                               folding.  [default: 1.0]
    @param pupil_plane_im      The GalSim.Image, NumPy array, or name of file containing the pupil
                               plane image, to be used instead of generating one based on the
                               obscuration and strut parameters.  [default: None]
    @param pupil_angle         If `pupil_plane_im` is not None, rotation angle for the pupil plane
                               (positive in the counter-clockwise direction).  Must be an Angle
                               instance. [default: 0. * galsim.degrees]
    @param pupil_plane_scale   Sampling interval in meters to use for the pupil plane array.  In
                               most cases, it's a good idea to leave this as None, in which case
                               GalSim will attempt to find a good value automatically.  The
                               exception is when specifying the pupil arrangement via an image, in
                               which case this keyword can be used to indicate the sampling of that
                               image.  See also `pad_factor` for adjusting the pupil sampling scale.
                               [default: None]
    @param pupil_plane_size    Size in meters to use for the pupil plane array.  In most cases, it's
                               a good idea to leave this as None, in which case GalSim will attempt
                               to find a good value automatically.  See also `oversampling` for
                               adjusting the pupil size.  [default: None]
    """
    def __init__(self, screen_list, lam, t0=0.0, exptime=0.0, time_step=0.025, flux=1.0, aper=None,
                 theta=(0.0*galsim.arcmin, 0.0*galsim.arcmin), interpolant=None,
                 scale_unit=galsim.arcsec, ii_pad_factor=4., suppress_warning=False,
                 geometric_shooting=True, gsparams=None,
                 _bar=None, _force_stepk=None, _force_maxk=None, **kwargs):
        # Hidden `_bar` kwarg can be used with astropy.console.utils.ProgressBar to print out a
        # progress bar during long calculations.

        self._screen_list = screen_list
        self.t0 = float(t0)
        self.lam = float(lam)
        self.exptime = float(exptime)
        self.time_step = float(time_step)
        if aper is None:
            # Check here for diameter.
            if 'diam' not in kwargs:
                raise ValueError("Diameter required if aperture not specified directly.")
            aper = Aperture(lam=lam, screen_list=self._screen_list, gsparams=gsparams, **kwargs)
        self.aper = aper
        if not isinstance(theta[0], galsim.Angle) or not isinstance(theta[1], galsim.Angle):
            raise TypeError("theta must be 2-tuple of galsim.Angle's.")
        self.theta = theta
        self.interpolant = interpolant
        if isinstance(scale_unit, str):
            scale_unit = galsim.angle.get_angle_unit(scale_unit)
        self.scale_unit = scale_unit
        self._gsparams = gsparams
        self.scale = aper._sky_scale(self.lam, self.scale_unit)

        self._serialize_stepk = _force_stepk
        self._serialize_maxk = _force_maxk

        # Difference between serialize_maxk and force_maxk in InterpolatedImage is a factor of
        # scale.
        if self._serialize_stepk is not None:
            self._serialize_stepk *= self.scale
        if self._serialize_maxk is not None:
            self._serialize_maxk *= self.scale

        self.img = np.zeros(self.aper.illuminated.shape, dtype=np.float64)

        if self.exptime < 0:
            raise ValueError("Cannot integrate PSF for negative time.")

        self._ii_pad_factor = ii_pad_factor

        self._bar = _bar
        self._flux = flux
        self._suppress_warning = suppress_warning
        self._geometric_shooting = geometric_shooting

        # Need to put in a placeholder SBProfile so that calls to, for example,
        # self.stepk, still work.
        array = np.array([[self._flux]], dtype=np.float)
        bounds = galsim._BoundsI(1, 1, 1, 1)
        wcs = galsim.PixelScale(self.scale)
        image = galsim._Image(array, bounds, wcs)
        dummy_interpolant = 'delta' # so wavefront gradient photon-shooting works.
        self._dummy_obj = galsim.InterpolatedImage(
                image, pad_factor=1.0, x_interpolant=dummy_interpolant,
                _serialize_stepk=self._serialize_stepk,
                _serialize_maxk=self._serialize_maxk)
        self._sbp = self._dummy_obj._sbp

        self._screen_list._delayCalculation(self)

    @property
    def flux(self):
        return self._flux

    def __str__(self):
        return ("galsim.PhaseScreenPSF(%s, lam=%s, exptime=%s)" %
                (self._screen_list, self.lam, self.exptime))

    def __repr__(self):
        outstr = ("galsim.PhaseScreenPSF(%r, lam=%r, exptime=%r, flux=%r, aper=%r, theta=%r, " +
                  "interpolant=%r, scale_unit=%r, gsparams=%r)")
        return outstr % (self._screen_list, self.lam, self.exptime, self.flux, self.aper,
                         self.theta, self.interpolant, self.scale_unit, self.gsparams)

    def __eq__(self, other):
        # Even if two PSFs were generated with different sets of parameters, they will act
        # identically if their img, interpolant, stepk, maxk, pad_factor, and gsparams match.
        return (isinstance(other, PhaseScreenPSF) and
                self._screen_list == other._screen_list and
                self.lam == other.lam and
                self.aper == other.aper and
                self.t0 == other.t0 and
                self.exptime == other.exptime and
                self.time_step == other.time_step and
                self._flux == other._flux and
                self.interpolant == other.interpolant and
                self._serialize_stepk == other._serialize_stepk and
                self._serialize_maxk == other._serialize_maxk and
                self._ii_pad_factor == other._ii_pad_factor and
                self.gsparams == other.gsparams)

    def __hash__(self):
        return hash(("galsim.PhaseScreenPSF", tuple(self._screen_list), self.lam, self.aper,
                     self.t0, self.exptime, self.time_step, self._flux, self.interpolant,
                     self._serialize_stepk, self._serialize_maxk, self._ii_pad_factor,
                     self.gsparams))

    def _prepareDraw(self):
        # Trigger delayed computation of all pending PSFs.
        self._screen_list._prepareDraw()

    # A few items which need the InterpolatedImage to have been prepared before accessing.
    @property
    def maxk(self):
        """The value of k beyond which aliasing can be neglected.
        """
        self._prepareDraw()
        return self.ii.maxk

    @property
    def nyquist_scale(self):
        """The Image pixel spacing that does not alias maxk.
        """
        # Use this instead of self.ii.nyquistScale() so we don't need to _prepareDraw when
        # photon-shooting into an automatically-sized image.
        return np.pi/self.aper._maxK(self.lam, self.scale_unit)

    @property
    def stepk(self):
        """The sampling in k space necessary to avoid folding of image in x space.
        """
        self._prepareDraw()
        return self.ii.stepk

    @property
    def centroid(self):
        """The (x, y) centroid of an object as a Position.
        """
        self._prepareDraw()
        return self.ii.centroid

    @property
    def max_sb(self):
        """An estimate of the maximum surface brightness of the object.

        Some profiles will return the exact peak SB, typically equal to the value of
        obj.xValue(obj.centroid).  However, not all profiles (e.g. Convolution) know how to
        calculate this value without just drawing the image and checking what the maximum value is.
        Clearly, this would be inefficient, so in these cases, some kind of estimate is returned,
        which will generally be conservative on the high side.

        This routine is mainly used by the photon shooting process, where an overestimate of
        the maximum surface brightness is acceptable.

        Note, for negative-flux profiles, this will return the absolute value of the most negative
        surface brightness.  Technically, it is an estimate of the maximum deviation from zero,
        rather than the maximum value.  For most profiles, these are the same thing.
        """
        self._prepareDraw()
        return self.ii.max_sb

    def _step(self):
        """Compute the current instantaneous PSF and add it to the developing integrated PSF."""
        u = self.aper.u[self.aper.illuminated]
        v = self.aper.v[self.aper.illuminated]
        wf = self._screen_list._wavefront(u, v, None, self.theta)
        expwf = np.exp((2j*np.pi/self.lam) * wf)
        expwf_grid = np.zeros_like(self.aper.illuminated, dtype=np.complex128)
        expwf_grid[self.aper.illuminated] = expwf
        ftexpwf = galsim.fft.fft2(expwf_grid, shift_in=True, shift_out=True)
        self.img += np.abs(ftexpwf)**2

    def _finalize(self):
        """Take accumulated integrated PSF image and turn it into a proper GSObject."""
        self.img *= self._flux / self.img.sum(dtype=float)
        b = galsim._BoundsI(1,self.aper.npix,1,self.aper.npix)
        self.img = galsim._Image(self.img, b, galsim.PixelScale(self.scale))

        self.ii = galsim.InterpolatedImage(
                self.img, x_interpolant=self.interpolant,
                _serialize_stepk=self._serialize_stepk, _serialize_maxk=self._serialize_maxk,
                pad_factor=self._ii_pad_factor,
                use_true_center=False, gsparams=self._gsparams)

        self._sbp = self.ii._sbp

        if not self._suppress_warning:
            specified_stepk = 2*np.pi/(self.img.array.shape[0]*self.scale)
            observed_stepk = self.ii.stepk

            if observed_stepk < specified_stepk:
                import warnings
                warnings.warn(
                    "The calculated stepk (%g) for PhaseScreenPSF is smaller "%observed_stepk +
                    "than what was used to build the wavefront (%g). "%specified_stepk +
                    "This could lead to aliasing problems. " +
                    "Increasing pad_factor is recommended.")

    def __getstate__(self):
        # Finish calculating before pickling.
        self._prepareDraw()
        d = self.__dict__.copy()
        # The SBProfile is picklable, but it is pretty inefficient, due to the large images being
        # written as a string.  Better to pickle the image and remake the InterpolatedImage.
        del d['_sbp']
        del d['ii']
        d.pop('_dummy_obj',None)
        return d

    def __setstate__(self, d):
        self.__dict__ = d
        self.ii = galsim.InterpolatedImage(self.img, x_interpolant=self.interpolant,
                                           use_true_center=False,
                                           pad_factor=self._ii_pad_factor,
                                           _serialize_stepk=self._serialize_stepk,
                                           _serialize_maxk=self._serialize_maxk,
                                           gsparams=self._gsparams)
        self._sbp = self.ii._sbp

    def shoot(self, n_photons, rng=None):
        """Shoot photons into a PhotonArray.

        @param n_photons    The number of photons to use for photon shooting.
        @param rng          If provided, a random number generator to use for photon shooting,
                            which may be any kind of BaseDeviate object.  If `rng` is None, one
                            will be automatically created, using the time as a seed.
                            [default: None]
        @returns PhotonArray.
        """
        if not self._geometric_shooting:
            self._prepareDraw()
            return self.ii.shoot(n_photons, rng)

        ud = galsim.UniformDeviate(rng)

        t = np.empty((n_photons,), dtype=float)
        ud.generate(t)
        t *= self.exptime
        t += self.t0
        u = self.aper.u[self.aper.illuminated]
        v = self.aper.v[self.aper.illuminated]
        pick = np.empty((n_photons,), dtype=float)
        ud.generate(pick)
        pick *= len(u)
        pick = pick.astype(int)
        u = u[pick]
        v = v[pick]

        x, y = self._screen_list._wavefront_gradient(u, v, t, self.theta)
        x *= 1e-9 * 206265  # convert wavefront gradient from nm/m to arcsec.
        y *= 1e-9 * 206265

        photon_array = galsim._galsim.PhotonArray(n_photons)
        photon_array.x = x
        photon_array.y = y
        photon_array.flux = self._flux/n_photons
        return photon_array


class OpticalPSF(GSObject):
    """A class describing aberrated PSFs due to telescope optics.  Its underlying implementation
    uses an InterpolatedImage to characterize the profile.

    The diffraction effects are characterized by the diffraction angle, which is a function of the
    ratio lambda / D, where lambda is the wavelength of the light and D is the diameter of the
    telescope.  The natural unit for this value is radians, which is not normally a convenient
    unit to use for other GSObject dimensions.  Assuming that the other sky coordinates you are
    using are all in arcsec (e.g. the pixel scale when you draw the image, the size of the galaxy,
    etc.), then you should convert this to arcsec as well:

        >>> lam = 700  # nm
        >>> diam = 4.0    # meters
        >>> lam_over_diam = (lam * 1.e-9) / diam  # radians
        >>> lam_over_diam *= 206265  # Convert to arcsec
        >>> psf = galsim.OpticalPSF(lam_over_diam, ...)

    To make this process a bit simpler, we recommend instead providing the wavelength and diameter
    separately using the parameters `lam` (in nm) and `diam` (in m).  GalSim will then convert this
    to any of the normal kinds of angular units using the `scale_unit` parameter:

        >>> psf = galsim.OpticalPSF(lam=lam, diam=diam, scale_unit=galsim.arcsec, ...)

    When drawing images, the scale_unit should match the unit used for the pixel scale or the WCS.
    e.g. in this case, a pixel scale of 0.2 arcsec/pixel would be specified as `pixel_scale=0.2`.

    Input aberration coefficients are assumed to be supplied in units of wavelength, and correspond
    to the Zernike polynomials in the Noll convention defined in
    Noll, J. Opt. Soc. Am. 66, 207-211(1976).  For a brief summary of the polynomials, refer to
    http://en.wikipedia.org/wiki/Zernike_polynomials#Zernike_polynomials.  By default, the
    aberration coefficients indicate the amplitudes of _circular_ Zernike polynomials, which are
    orthogonal over a circle.  If you would like the aberration coefficients to instead be
    interpretted as the amplitudes of _annular_ Zernike polynomials, which are orthogonal over an
    annulus (see Mahajan, J. Opt. Soc. Am. 71, 1 (1981)), set the `annular_zernike` keyword argument
    to True.

    There are two ways to specify the geometry of the pupil plane, i.e., the obscuration disk size
    and the areas that will be illuminated outside of it.  The first way is to use keywords that
    specify the size of the obscuration, and the nature of the support struts holding up the
    secondary mirror (or prime focus cage, etc.).  These are taken to be rectangular obscurations
    extending from the outer edge of the pupil to the outer edge of the obscuration disk (or the
    pupil center if `obscuration = 0.`).  You can specify how many struts there are (evenly spaced
    in angle), how thick they are as a fraction of the pupil diameter, and what angle they start at
    relative to the positive y direction.

    The second way to specify the pupil plane configuration is by passing in an image of it.  This
    can be useful for example if the struts are not evenly spaced or are not radially directed, as
    is assumed by the simple model for struts described above.  In this case, keywords related to
    struts are ignored; moreover, the `obscuration` keyword is used to ensure that the images are
    properly sampled (so it is still needed), but the keyword is then ignored when using the
    supplied image of the pupil plane.  Note that for complicated pupil configurations, it may be
    desireable to increase `pad_factor` for more fidelity at the expense of slower running time.
    The `pupil_plane_im` that is passed in can be rotated during internal calculations by specifying
    a `pupil_angle` keyword.

    If you choose to pass in a pupil plane image, it must be a square array in which the image of
    the pupil is centered.  The areas that are illuminated should have some value >0, and the other
    areas should have a value of precisely zero.  Based on what the OpticalPSF class thinks is the
    required sampling to make the PSF image, the image that is passed in of the pupil plane might be
    zero-padded during internal calculations.  The pixel scale of the pupil plane can be specified
    in one of three ways.  In descending order of priority, these are:
      1.  The `pupil_plane_scale` keyword argument (units are meters).
      2.  The `pupil_plane_im.scale` attribute (units are meters).
      3.  If (1) and (2) are both None, then the scale will be inferred by assuming that the
          illuminated pixel farthest from the image center is at a physical distance of self.diam/2.
    Note that if the scale is specified by either (1) or (2) above (which always includes specifying
    the pupil_plane_im as a filename, since the default scale then will be 1.0), then the
    lam_over_diam keyword must not be used, but rather the lam and diam keywords are required
    separately.  Finally, to ensure accuracy of calculations using a pupil plane image, we recommend
    sampling it as finely as possible.

    Initialization
    --------------

    As described above, either specify the lam/diam ratio directly in arbitrary units:

        >>> optical_psf = galsim.OpticalPSF(lam_over_diam=lam_over_diam, defocus=0., ...)

    or, use separate keywords for the telescope diameter and wavelength in meters and nanometers,
    respectively:

        >>> optical_psf = galsim.OpticalPSF(lam=lam, diam=diam, defocus=0., ...)

    Either of these options initializes `optical_psf` as an OpticalPSF instance.

    @param lam_over_diam    Lambda / telescope diameter in the physical units adopted for `scale`
                            (user responsible for consistency).  Either `lam_over_diam`, or `lam`
                            and `diam`, must be supplied.
    @param lam              Lambda (wavelength) in units of nanometers.  Must be supplied with
                            `diam`, and in this case, image scales (`scale`) should be specified in
                            units of `scale_unit`.
    @param diam             Telescope diameter in units of meters.  Must be supplied with
                            `lam`, and in this case, image scales (`scale`) should be specified in
                            units of `scale_unit`.
    @param tip              Tip in units of incident light wavelength. [default: 0]
    @param tilt             Tilt in units of incident light wavelength. [default: 0]
    @param defocus          Defocus in units of incident light wavelength. [default: 0]
    @param astig1           Astigmatism (like e2) in units of incident light wavelength.
                            [default: 0]
    @param astig2           Astigmatism (like e1) in units of incident light wavelength.
                            [default: 0]
    @param coma1            Coma along y in units of incident light wavelength. [default: 0]
    @param coma2            Coma along x in units of incident light wavelength. [default: 0]
    @param trefoil1         Trefoil (one of the arrows along y) in units of incident light
                            wavelength. [default: 0]
    @param trefoil2         Trefoil (one of the arrows along x) in units of incident light
                            wavelength. [default: 0]
    @param spher            Spherical aberration in units of incident light wavelength.
                            [default: 0]
    @param aberrations      Optional keyword, to pass in a list, tuple, or NumPy array of
                            aberrations in units of reference wavelength (ordered according to
                            the Noll convention), rather than passing in individual values for each
                            individual aberration.  Note that aberrations[1] is piston (and not
                            aberrations[0], which is unused.)  This list can be arbitrarily long to
                            handle Zernike polynomial aberrations of arbitrary order.
    @param annular_zernike  Boolean indicating that aberrations specify the amplitudes of annular
                            Zernike polynomials instead of circular Zernike polynomials.
                            [default: False]
    @param aper             Aperture object to use when creating PSF.  [default: None]
    @param circular_pupil   Adopt a circular pupil?  [default: True]
    @param obscuration      Linear dimension of central obscuration as fraction of pupil linear
                            dimension, [0., 1.). This should be specified even if you are providing
                            a `pupil_plane_im`, since we need an initial value of obscuration to use
                            to figure out the necessary image sampling. [default: 0]
    @param interpolant      Either an Interpolant instance or a string indicating which interpolant
                            should be used.  Options are 'nearest', 'sinc', 'linear', 'cubic',
                            'quintic', or 'lanczosN' where N should be the integer order to use.
                            [default: galsim.Quintic()]
    @param oversampling     Optional oversampling factor for the InterpolatedImage. Setting
                            `oversampling < 1` will produce aliasing in the PSF (not good).
                            Usually `oversampling` should be somewhat larger than 1.  1.5 is
                            usually a safe choice.  [default: 1.5]
    @param pad_factor       Additional multiple by which to zero-pad the PSF image to avoid folding
                            compared to what would be employed for a simple Airy.  Note that
                            `pad_factor` may need to be increased for stronger aberrations, i.e.
                            those larger than order unity.  [default: 1.5]
    @param ii_pad_factor    Zero-padding factor by which to extend the image of the PSF when
                            creating the `InterpolatedImage`.  See the `InterpolatedImage` docstring
                            for more details.  [default: 4.]
    @param suppress_warning If `pad_factor` is too small, the code will emit a warning telling you
                            its best guess about how high you might want to raise it.  However,
                            you can suppress this warning by using `suppress_warning=True`.
                            [default: False]
    @param geometric_shooting  If True, then when drawing using photon shooting, use geometric
                            optics approximation where the photon angles are derived from the
                            phase screen gradient.  If False, then first draw using Fourier
                            optics and then shoot from the derived InterpolatedImage.
                            [default: False]
    @param flux             Total flux of the profile. [default: 1.]
    @param nstruts          Number of radial support struts to add to the central obscuration.
                            [default: 0]
    @param strut_thick      Thickness of support struts as a fraction of pupil diameter.
                            [default: 0.05]
    @param strut_angle      Angle made between the vertical and the strut starting closest to it,
                            defined to be positive in the counter-clockwise direction; must be an
                            Angle instance. [default: 0. * galsim.degrees]
    @param pupil_plane_im   The GalSim.Image, NumPy array, or name of file containing the pupil
                            plane image, to be used instead of generating one based on the
                            obscuration and strut parameters.  [default: None]
    @param pupil_angle      If `pupil_plane_im` is not None, rotation angle for the pupil plane
                            (positive in the counter-clockwise direction).  Must be an Angle
                            instance. [default: 0. * galsim.degrees]
    @param pupil_plane_scale Sampling interval in meters to use for the pupil plane array.  In
                            most cases, it's a good idea to leave this as None, in which case
                            GalSim will attempt to find a good value automatically.  The
                            exception is when specifying the pupil arrangement via an image, in
                            which case this keyword can be used to indicate the sampling of that
                            image.  See also `pad_factor` for adjusting the pupil sampling scale.
                            [default: None]
    @param pupil_plane_size Size in meters to use for the pupil plane array.  In most cases, it's
                            a good idea to leave this as None, in which case GalSim will attempt
                            to find a good value automatically.  See also `oversampling` for
                            adjusting the pupil size.  [default: None]
    @param scale_unit       Units to use for the sky coordinates when calculating lam/diam if these
                            are supplied separately.  Should be either a galsim.AngleUnit or a
                            string that can be used to construct one (e.g., 'arcsec', 'radians',
                            etc.).  [default: galsim.arcsec]
    @param gsparams         An optional GSParams argument.  See the docstring for GSParams for
                            details. [default: None]

    Methods
    -------

    There are no additional methods for OpticalPSF beyond the usual GSObject methods.
    """
    _req_params = {}
    _opt_params = {
        "diam": float,
        "defocus": float,
        "astig1": float,
        "astig2": float,
        "coma1": float,
        "coma2": float,
        "trefoil1": float,
        "trefoil2": float,
        "spher": float,
        "annular_zernike": bool,
        "circular_pupil": bool,
        "obscuration": float,
        "oversampling": float,
        "pad_factor": float,
        "suppress_warning": bool,
        "max_size": float,
        "interpolant": str,
        "flux": float,
        "nstruts": int,
        "strut_thick": float,
        "strut_angle": galsim.Angle,
        "pupil_plane_im": str,
        "pupil_angle": galsim.Angle,
        "pupil_plane_scale": float,
        "pupil_plane_size": float,
        "scale_unit": str}
    _single_params = [{"lam_over_diam": float, "lam": float}]
    _takes_rng = False

    def __init__(self, lam_over_diam=None, lam=None, diam=None, tip=0., tilt=0., defocus=0.,
                 astig1=0., astig2=0., coma1=0., coma2=0., trefoil1=0., trefoil2=0., spher=0.,
                 aberrations=None, annular_zernike=False,
                 aper=None, circular_pupil=True, obscuration=0., interpolant=None,
                 oversampling=1.5, pad_factor=1.5, ii_pad_factor=4., flux=1.,
                 nstruts=0, strut_thick=0.05, strut_angle=0.*galsim.degrees,
                 pupil_plane_im=None, pupil_plane_scale=None, pupil_plane_size=None,
                 pupil_angle=0.*galsim.degrees, scale_unit=galsim.arcsec, gsparams=None,
                 _force_maxk=None, _force_stepk=None,
                 suppress_warning=False, geometric_shooting=False, max_size=None):
        if max_size is not None: # pragma: no cover
            from .deprecated import depr
            depr('max_size', 1.4, '',
                 "The max_size keyword has been removed.  In its place, the pad_factor keyword can"
                 "be used to adjust the size of the internal InterpolatedImage.")

        if isinstance(scale_unit, str):
            scale_unit = galsim.angle.get_angle_unit(scale_unit)
        # Need to handle lam/diam vs. lam_over_diam here since lam by itself is needed for
        # OpticalScreen.
        if lam_over_diam is not None:
            if lam is not None or diam is not None:
                raise TypeError("If specifying lam_over_diam, then do not specify lam or diam")
            # For combination of lam_over_diam and pupil_plane_im with a specified scale, it's
            # tricky to determine the actual diameter of the pupil needed by Aperture.  So for now,
            # we just disallow this combination.  Please feel free to raise an issue at
            # https://github.com/GalSim-developers/GalSim/issues if you need this functionality.
            if pupil_plane_im is not None:
                if isinstance(pupil_plane_im, str):
                    # Filename, therefore specific scale exists.
                    raise TypeError("If specifying lam_over_diam, then do not "
                                    "specify pupil_plane_im as a filename.")
                elif (isinstance(pupil_plane_im, galsim.Image)
                      and pupil_plane_im.scale is not None):
                    raise TypeError("If specifying lam_over_diam, then do not specify "
                                    "pupil_plane_im with definite scale attribute.")
                elif pupil_plane_scale is not None:
                    raise TypeError("If specifying lam_over_diam, then do not specify "
                                    "pupil_plane_scale.")
            lam = 500.  # Arbitrary
            diam = lam*1.e-9 / lam_over_diam * galsim.radians / scale_unit
        else:
            if lam is None or diam is None:
                raise TypeError("If not specifying lam_over_diam, then specify lam AND diam")

        # Make the optical screen.
        optics_screen = galsim.OpticalScreen(
                diam=diam, defocus=defocus, astig1=astig1, astig2=astig2, coma1=coma1, coma2=coma2,
                trefoil1=trefoil1, trefoil2=trefoil2, spher=spher, aberrations=aberrations,
                obscuration=obscuration, annular_zernike=annular_zernike, lam_0=lam)
        self._screens = galsim.PhaseScreenList(optics_screen)

        # Make the aperture.
        if aper is None:
            aper = galsim.Aperture(
                    diam, lam=lam, circular_pupil=circular_pupil, obscuration=obscuration,
                    nstruts=nstruts, strut_thick=strut_thick, strut_angle=strut_angle,
                    oversampling=oversampling, pad_factor=pad_factor,
                    pupil_plane_im=pupil_plane_im, pupil_angle=pupil_angle,
                    pupil_plane_scale=pupil_plane_scale, pupil_plane_size=pupil_plane_size,
                    gsparams=gsparams)
            self.obscuration = obscuration
        else:
            if hasattr(aper, '_obscuration'):
                self.obscuration = aper._obscuration
            else:
                self.obscuration = 0.0

        # Save for pickling
        self._lam = lam
        self._flux = flux
        self._interpolant = interpolant
        self._scale_unit = scale_unit
        self._gsparams = gsparams
        self._suppress_warning = suppress_warning
        self._geometric_shooting = geometric_shooting
        self._aper = aper
        self._force_maxk = _force_maxk
        self._force_stepk = _force_stepk
        self._ii_pad_factor = ii_pad_factor

        # Finally, put together to make the PSF.
        self._psf = galsim.PhaseScreenPSF(self._screens, lam=self._lam, flux=self._flux,
                                          aper=aper, interpolant=self._interpolant,
                                          scale_unit=self._scale_unit, gsparams=self._gsparams,
                                          suppress_warning=self._suppress_warning,
                                          geometric_shooting=self._geometric_shooting,
                                          _force_maxk=_force_maxk, _force_stepk=_force_stepk,
                                          ii_pad_factor=ii_pad_factor)

        self._psf._prepareDraw()  # No need to delay an OpticalPSF.
        self._sbp = self._psf._sbp

    @property
    def flux(self):
        return self._flux

    def __str__(self):
        screen = self._psf._screen_list[0]
        s = "galsim.OpticalPSF(lam=%s, diam=%s" % (screen.lam_0, self._psf.aper.diam)
        if any(screen.aberrations):
            s += ", aberrations=[" + ",".join(str(ab) for ab in screen.aberrations) + "]"
        if hasattr(self._psf.aper, '_circular_pupil'):
            s += self._psf.aper._geometry_str()
        if screen.annular_zernike:
            s += ", annular_zernike=True"
            s += ", obscuration=%r"%self.obscuration
        if self._flux != 1.0:
            s += ", flux=%s" % self._flux
        s += ")"
        return s

    def __repr__(self):
        screen = self._psf._screen_list[0]
        s = "galsim.OpticalPSF(lam=%r, diam=%r" % (self._lam, self._psf.aper.diam)
        s += ", aper=%r"%self._psf.aper
        if any(screen.aberrations):
            s += ", aberrations=[" + ",".join(repr(ab) for ab in screen.aberrations) + "]"
        if screen.annular_zernike:
            s += ", annular_zernike=True"
            s += ", obscuration=%r"%self.obscuration
        if self._flux != 1.0:
            s += ", flux=%r" % self._flux
        if self._force_maxk is not None:
            s += ", _force_maxk=%r" % self._force_maxk
        if self._force_stepk is not None:
            s += ", _force_stepk=%r" % self._force_stepk
        if self._ii_pad_factor is not None:
            s += ", ii_pad_factor=%r" % self._ii_pad_factor
        s += ")"
        return s

    def __eq__(self, other):
        return (isinstance(other, galsim.OpticalPSF) and
                self._lam == other._lam and
                self._aper == other._aper and
                self._psf._screen_list == other._psf._screen_list and
                self._flux == other._flux and
                self._interpolant == other._interpolant and
                self._scale_unit == other._scale_unit and
                self._force_stepk == other._force_stepk and
                self._force_maxk == other._force_maxk and
                self._ii_pad_factor == other._ii_pad_factor and
                self._gsparams == other._gsparams)

    def __hash__(self):
        return hash(("galsim.OpticalPSF", self._lam, self._aper, self._psf._screen_list[0],
                     self._flux, self._interpolant, self._scale_unit, self._force_stepk,
                     self._force_maxk, self._ii_pad_factor, self._gsparams))

    def __getstate__(self):
        # The SBProfile is picklable, but it is pretty inefficient, due to the large images being
        # written as a string.  Better to pickle the psf and remake the PhaseScreenPSF.
        d = self.__dict__.copy()
        d['aper'] = d['_psf'].aper
        del d['_sbp']
        del d['_psf']
        return d

    def __setstate__(self, d):
        self.__dict__ = d
        aper = self.__dict__.pop('aper')
        self._psf = galsim.PhaseScreenPSF(self._screens, lam=self._lam, flux=self._flux,
                                          aper=aper, interpolant=self._interpolant,
                                          scale_unit=self._scale_unit, gsparams=self._gsparams,
                                          suppress_warning=self._suppress_warning,
                                          _force_maxk=self._force_maxk,
                                          _force_stepk=self._force_stepk,
                                          ii_pad_factor=self._ii_pad_factor)
        self._psf._prepareDraw()
        self._sbp = self._psf._sbp

    def shoot(self, n_photons, rng=None):
        """Shoot photons into a PhotonArray.

        @param n_photons    The number of photons to use for photon shooting.
        @param rng          If provided, a random number generator to use for photon shooting,
                            which may be any kind of BaseDeviate object.  If `rng` is None, one
                            will be automatically created, using the time as a seed.
                            [default: None]
        @returns PhotonArray.
        """
        return self._psf.shoot(n_photons, rng)
