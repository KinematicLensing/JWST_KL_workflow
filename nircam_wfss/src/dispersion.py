"""
dispersion.py
-------------
Grism dispersion and spectral-trace polynomial models, plus the
``grism_conf_preparation`` helper that converts a direct-image pixel
position to arrays of (dx, dy, wavelength) for spectral extraction.

Polynomial model notes
~~~~~~~~~~~~~~~~~~~~~~
These are custom (non-grismconf) polynomial fits to the NIRCam grism
dispersion relation, calibrated by F. Sun.  Two functions are defined:

* ``fit_disp_order32`` – maps (x_pix, y_pix, wavelength) → dx (pixel offset
  along the dispersion axis). Third-order in wavelength, second-order in
  pixel position.

* ``fit_disp_order23`` – maps (x_pix, y_pix, dx) → dy (pixel offset in the
  cross-dispersion direction).  Second-order in dx, third-order in pixel
  position.

Their parameter vectors are stored in calibration files and read by the
pipeline before calling ``grism_conf_preparation``.
"""

from __future__ import annotations

import numpy as np
from scipy import interpolate
import os
from astropy.io import ascii


# ---------------------------------------------------------------------------
# Dispersion polynomial: dx = f(x_pix, y_pix, lambda)
# ---------------------------------------------------------------------------

def fit_disp_order32(
    data: np.ndarray,
    a01: float, a02: float, a03: float, a04: float, a05: float, a06: float,
    b01: float, b02: float, b03: float, b04: float, b05: float, b06: float,
    c01: float, c02: float, c03: float,
    d01: float,
) -> np.ndarray:
    """
    Dispersion model: pixel offset dx as a function of position and wavelength.

    Parameterisation: third-order in (wavelength − 3.95 µm),
    second-order in (x_pix − 1024) and (y_pix − 1024).

    Parameters
    ----------
    data:
        Shape (3, N) array where
        ``data[0]`` = x pixel,
        ``data[1]`` = y pixel,
        ``data[2]`` = wavelength (µm).
    a01 … d01:
        Polynomial coefficients (30 total, matching the calibration file
        convention from F. Sun).

    Returns
    -------
    dx:
        Pixel offset along the dispersion axis.
    """
    xpix, ypix, dx = data[0] - 1024, data[1] - 1024, data[2] - 3.95
    return (
        (
            a01
            + (a02 * xpix + a03 * ypix)
            + (a04 * xpix**2 + a05 * xpix * ypix + a06 * ypix**2)
        )
        + (
            b01
            + (b02 * xpix + b03 * ypix)
            + (b04 * xpix**2 + b05 * xpix * ypix + b06 * ypix**2)
        ) * dx
        + (
            c01
            + (c02 * xpix + c03 * ypix)
        ) * dx**2
        + d01 * dx**3
    )


# ---------------------------------------------------------------------------
# Trace polynomial: dy = f(x_pix, y_pix, dx)
# ---------------------------------------------------------------------------

def fit_disp_order23(
    data: np.ndarray,
    a01: float, a02: float, a03: float, a04: float, a05: float,
    a06: float, a07: float, a08: float, a09: float, a10: float,
    b01: float, b02: float, b03: float, b04: float, b05: float,
    b06: float, b07: float, b08: float, b09: float, b10: float,
    c01: float, c02: float, c03: float, c04: float, c05: float,
    c06: float, c07: float, c08: float, c09: float, c10: float,
) -> np.ndarray:
    """
    Trace model: cross-dispersion pixel offset dy as a function of position
    and dispersion offset dx.

    Parameterisation: second-order in dx, third-order in
    (x_pix − 1024) and (y_pix − 1024).

    Parameters
    ----------
    data:
        Shape (3, N) array where
        ``data[0]`` = x pixel,
        ``data[1]`` = y pixel,
        ``data[2]`` = dx (dispersion offset in pixels).
    a01 … c10:
        30 polynomial coefficients.

    Returns
    -------
    dy:
        Pixel offset in the cross-dispersion direction.
    """
    xpix, ypix, dx = data[0] - 1024, data[1] - 1024, data[2]
    poly3 = lambda p, q, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10: (  # noqa: E731
        c1
        + (c2 * p + c3 * q)
        + (c4 * p**2 + c5 * p * q + c6 * q**2)
        + (c7 * p**3 + c8 * p**2 * q + c9 * p * q**2 + c10 * q**3)
    )
    return (
        poly3(xpix, ypix, a01, a02, a03, a04, a05, a06, a07, a08, a09, a10)
        + poly3(xpix, ypix, b01, b02, b03, b04, b05, b06, b07, b08, b09, b10) * dx
        + poly3(xpix, ypix, c01, c02, c03, c04, c05, c06, c07, c08, c09, c10) * dx**2
    )


# ---------------------------------------------------------------------------
# High-level helper: prepare dxs, dys, wavelengths for one source position
# ---------------------------------------------------------------------------

def grism_conf_preparation(
    x0: float = 1024.0,
    y0: float = 1024.0,
    pupil: str = "R",
    fit_opt_fit: np.ndarray | None = None,
    w_opt: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute spectral trace offsets and wavelengths for a given source position.

    Uses the custom polynomial calibration (``fit_disp_order32`` for the
    wavelength-to-pixel mapping and ``fit_disp_order23`` for the trace shape)
    rather than the ``grismconf`` library, giving finer control over the
    calibration coefficients.

    Parameters
    ----------
    x0, y0:
        Source position in the direct image (pixels).
    pupil:
        Grism pupil wheel position: ``'R'`` (GRISMR, disperses along X) or
        ``'C'`` (GRISMC, disperses along Y).
    fit_opt_fit:
        1-D array of 30 trace polynomial coefficients for ``fit_disp_order23``.
        Defaults to zeros (flat trace).
    w_opt:
        1-D array of 16 dispersion polynomial coefficients for
        ``fit_disp_order32``.  Defaults to zeros.

    Returns
    -------
    dxs:
        Pixel offsets along X from the direct-image position.
    dys:
        Pixel offsets along Y from the direct-image position.
    wavs:
        Wavelength array (µm) corresponding to each (dx, dy) step.
    """
    if fit_opt_fit is None:
        fit_opt_fit = np.zeros(30)
    if w_opt is None:
        w_opt = np.zeros(16)

    # Build wavelength → dx mapping over the NIRCam LW range
    wave_space = np.arange(2.39, 5.15, 0.01)
    disp_space = fit_disp_order32(
        np.vstack(
            (
                x0 * np.ones_like(wave_space),
                y0 * np.ones_like(wave_space),
                wave_space,
            )
        ),
        *w_opt,
    )
    # Invert: dx → wavelength
    sort_idx = np.argsort(disp_space)
    inverse_wave_disp = interpolate.UnivariateSpline(
        disp_space[sort_idx], wave_space[sort_idx], s=0, k=1
    )

    if pupil == "R":
        dxs = np.arange(int(np.min(disp_space)), int(np.max(disp_space))) - x0 % 1
        wavs = inverse_wave_disp(dxs)
        dys = fit_disp_order23(
            np.vstack(
                (x0 * np.ones_like(dxs), y0 * np.ones_like(dxs), dxs)
            ),
            *fit_opt_fit,
        )
    elif pupil == "C":
        dys = np.arange(int(np.min(disp_space)), int(np.max(disp_space))) - y0 % 1
        wavs = inverse_wave_disp(dys)
        dxs = fit_disp_order23(
            np.vstack(
                (x0 * np.ones_like(dys), y0 * np.ones_like(dys), dys)
            ),
            *fit_opt_fit,
        )
    else:
        raise ValueError("pupil must be 'R' or 'C'")

    return dxs, dys, wavs

class GrismConf:
    """
    Container for grism dispersion and trace polynomial coefficients and sensitivity curve.
    """
    __WRANGE__ = {
        'F444W': np.array([3.8, 5.1]),
        'F322W2': np.array([2.4, 4.1]),
        'F356W':  np.array([3.1, 4.0]),
        'F277W':  np.array([2.4, 3.1]),
    }
    __DISP_FILTER__ = {
        'F277W': 'F322W2',
        'F335M': 'F322W2',
        'F322W2': 'F322W2',
        'F356W': 'F322W2',
        'F360M': 'F322W2',
        'F410M': 'F444W',
        'F444W': 'F444W',
        'F480M': 'F444W',
    }

    def __init__(
        self,
        filter: str,
        config: dict,
    ) -> None:
        """ Initialize the GrismConf object with the given parameters.
        Parameters:
            filter:
                Grism filter name, e.g. 'F356W'.
            config:
                Configuration dictionary containing
        """
        self.filter = filter # filter for transmission
        # Default to full LW range if filter not recognized
        self.WRANGE = self.__WRANGE__.get(filter, np.array([2.4, 5.1]))
        self.disp_filter = self.__DISP_FILTER__.get(filter) # filter for dispersion
        self.config = config

        self.list_mod_pupil = []
        self.list_disp_coeff = []
        self.list_trace_coeff = []
        self.sensitivity = []
        for module in ['A', 'B']:
            for grism in ['R', 'C']:
                self.list_mod_pupil.append(f'{module}{grism}')
                self.list_trace_coeff.append(self._load_disp(module, grism)[0])
                self.list_disp_coeff.append(self._load_displ(module, grism)[0])
                self.sensitivity.append(self._load_sens(module, grism))

    def _load_disp(self, module, grism):
        """ Load the polynomial coefficients of spectral tracing model solution
            dy(x0, y0, dx), where dx = x_s - x0, dy = y_s - y0
        for the given module and grism. The coefficients are passed to `fit_disp_order23`.
        Parameters:
            module:
                NIRCam module ('A' or 'B').
            grism:
                Grism name ('R' or 'C').
        Returns:
            coeff_bestfit:
                1-D array of best-fitting trace polynomial coefficients for fit_disp_order23.
            coeff_error:
                1-D array of error of trace polynomial coefficients for fit_disp_order23.
        """
        path = os.path.join(
            self.config.cali_support_dir,
            "DISP_%s_mod%s_grism%s.dat" % (self.disp_filter, module, grism),
        )
        tb = ascii.read(path)
        return tb["col0"].data, tb["col1"].data
    
    def _load_displ(self, module, grism):
        """ Load the polynomial coefficients of dispersion wavelength solution
            dx(x0, y0, lambda_s)
        for the given module and grism. The coefficients are passed to `fit_disp_order32`.
        Parameters:
            module:
                NIRCam module ('A' or 'B').
            grism:
                Grism name ('R' or 'C').
        Returns:
            coeff_bestfit:
                1-D array of best-fitting trace polynomial coefficients for fit_disp_order32.
            coeff_error:
                1-D array of error of trace polynomial coefficients for fit_disp_order32.
        """
        path = os.path.join(
            self.config.cali_support_dir,
            "DISPL_mod%s_grism%s.dat" % (module, grism),
        )
        tb = ascii.read(path)
        return tb["col0"].data, tb["col1"].data
    
    def _load_sens(self, module, grism):
        """ Load the sensitivity curve for the given module and grism. 
        Parameters:
            module:
                NIRCam module ('A' or 'B').
            grism:
                Grism name ('R' or 'C').
        Returns:
            interpolate_sens:
                A callable function that takes wavelength as input and returns 
                the sensitivity (DN/s/Jy).
        """
        path = os.path.join(
            self.config.cali_support_dir,
            "%s_mod%s_grism%s_sensitivity.dat" % (self.filter, module, grism),
        )
        tb = ascii.read(path)
        interpolate_sens = interpolate.UnivariateSpline(tb["wavelength"], tb["DN/s/Jy"],
                                                        ext='zeros', k=1, s=1e2)
                                                        
        return interpolate_sens
    
    def get_disp_coeff(self, module, grism):
        """ Get the dispersion polynomial coefficients for the given module and grism.
        Parameters:
            module:
                NIRCam module ('A' or 'B').
            grism:
                Grism name ('R' or 'C').
        Returns:
            1-D array of dispersion polynomial coefficients for fit_disp_order32.
        """
        idx = self.list_mod_pupil.index(f'{module}{grism}')
        return self.list_disp_coeff[idx]

    def get_trace_coeff(self, module, grism):
        """ Get the trace polynomial coefficients for the given module and grism.
        Parameters:
            module:
                NIRCam module ('A' or 'B').
            grism:
                Grism name ('R' or 'C').
        Returns:
            1-D array of trace polynomial coefficients for fit_disp_order23.
        """
        idx = self.list_mod_pupil.index(f'{module}{grism}')
        return self.list_trace_coeff[idx]
    
    def get_sensitivity(self, module, grism):
        """ Get the sensitivity curve for the given module and grism.
        Parameters:
            module:
                NIRCam module ('A' or 'B').
            grism:
                Grism name ('R' or 'C').
        Returns:
            A callable function that takes wavelength as input and returns 
            the sensitivity (DN/s/Jy).
        """
        idx = self.list_mod_pupil.index(f'{module}{grism}')
        return self.sensitivity[idx]

# ---------------------------------------------------------------------------
# Convenience lambda functions used across the pipeline
# ---------------------------------------------------------------------------

#: Linear function  f(x) = k*x + b
linear = lambda x, k, b: x * k + b  # noqa: E731

#: Gaussian profile (integral = flux)
gauss = lambda x, x0, flux, fwhm: (  # noqa: E731
    (flux / fwhm / 1.064467)
    * np.exp(-((x - x0) ** 2) / (2 * (fwhm / 2.354820) ** 2))
)

#: Gaussian + linear continuum
gauss_cont_prof = lambda x, x0, flux, fwhm, k, b: (  # noqa: E731
    gauss(x, x0, flux, fwhm) + k * (x - 1024) + b
)
