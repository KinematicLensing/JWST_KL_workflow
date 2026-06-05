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
# High-level helper: get the detector position of a emission line, and its 
# local dispersion angle
# ---------------------------------------------------------------------------
def get_position_ang_dispang_at_wave(x0, y0, obswave, grism_conf, module, pupil, 
                                     delta_wave=0.005, wcs_transform=None, velosys=0.0):
    """
    Return the detector position of an emission line and its local dispersion angle.

    Parameters
    ----------
    x0, y0:
        Source position in the direct image (pixels).
    obswave:
        Observed wavelength of the emission line (µm).
    grism_conf:
        ``GrismConf`` instance with dispersion and trace polynomial coefficients.
    module:
        NIRCam module (``'A'`` or ``'B'``).
    pupil:
        Grism pupil: ``'R'`` (GRISMR, disperses along X) or ``'C'`` (GRISMC,
        disperses along Y).
    delta_wave:
        Small wavelength step (µm) used to compute the local dispersion angle
        via finite difference.
    wcs_transform:
        Optional ``[WCS_input, WCS_output]`` pair.  When provided, the
        dispersion angle is computed in the output WCS frame (e.g. a
        sky-aligned drizzled mosaic) rather than the native detector frame.
        Steps: predict line positions at ``obswave ± delta_wave`` in detector
        pixels, convert to sky via ``WCS_input``, project into the output
        frame via ``WCS_output``, then take ``arctan2(dy_out, dx_out)``.
    velosys:
        Systemic velocity (km/s) applied to the wavelength before solving the
        dispersion relation.  Default 0 (no shift).

    Returns
    -------
    xs, ys:
        Predicted detector position of the emission line (pixels), or
        ``None`` if ``obswave`` falls outside the modelled wavelength range.
    disp_ang:
        Local dispersion angle at the line position (radians).
    """
    ### Load dispersion solutions
    disp_coeff  = grism_conf.get_disp_coeff(module, pupil)
    trace_coeff = grism_conf.get_trace_coeff(module, pupil)
    dxs, dys, wavs = grism_conf_preparation(
        x0=x0, y0=y0, pupil=pupil,
        fit_opt_fit=trace_coeff, w_opt=disp_coeff,
    )
    wavs_bary = (1.0 + velosys / 299792458.0) * wavs

    if obswave < wavs_bary.min() or obswave > wavs_bary.max():
        return None

    interp_dx = interpolate.interp1d(
        wavs_bary, dxs, kind="linear", bounds_error=False, fill_value=np.nan
    )
    interp_dy = interpolate.interp1d(
        wavs_bary, dys, kind="linear", bounds_error=False, fill_value=np.nan
    )
    dx_line = float(interp_dx(obswave))
    dy_line = float(interp_dy(obswave))
    if np.isnan(dx_line) or np.isnan(dy_line):
        return None
    # predicted position of the emission line in the grism image
    xs = x0 + dx_line
    ys = y0 + dy_line

    ### Calculate the local dispersion angle at the grism WCS
    dx_lo = float(interp_dx(obswave - delta_wave))
    dy_lo = float(interp_dy(obswave - delta_wave))
    dx_hi = float(interp_dx(obswave + delta_wave))
    dy_hi = float(interp_dy(obswave + delta_wave))
    if wcs_transform:
        wcs_input, wcs_output = wcs_transform
        coords_lo = wcs_input.all_pix2world([[x0 + dx_lo, y0 + dy_lo]], 0)[0]
        coords_hi = wcs_input.all_pix2world([[x0 + dx_hi, y0 + dy_hi]], 0)[0]
        ox_lo, oy_lo = wcs_output.all_world2pix([[coords_lo[0], coords_lo[1]]], 0)[0]
        ox_hi, oy_hi = wcs_output.all_world2pix([[coords_hi[0], coords_hi[1]]], 0)[0]
        dispang = float(np.angle((ox_hi - ox_lo) + 1j * (oy_hi - oy_lo)))
    else:
        dispang = float(np.angle((dx_hi - dx_lo) + 1j * (dy_hi - dy_lo)))
    return xs, ys, dispang

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
    Container for grism dispersion, trace, and sensitivity calibrations.

    Loads calibration files for all four module/grism combinations (AR, AC, BR,
    BC) on construction.  Use ``get_disp_coeff``, ``get_trace_coeff``, and
    ``get_sensitivity`` to retrieve per-combination data.

    Parameters
    ----------
    filter:
        NIRCam filter name (e.g. ``'F444W'``, ``'F356W'``).  Determines the
        wavelength range and which sensitivity curve files are loaded.
    config:
        ``PipelineConfig`` instance; must expose ``cali_support_dir`` pointing
        to the directory that holds the calibration data files.

    Attributes
    ----------
    WRANGE:
        2-element array ``[wave_min, wave_max]`` (µm) for this filter.
    disp_filter:
        Filter name used to look up the dispersion coefficient files.  Some
        filters share coefficients (e.g. ``'F356W'`` → ``'F322W2'``).
    list_mod_pupil:
        List of module/grism keys: ``['AR', 'AC', 'BR', 'BC']``.
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
        """Initialize all calibration arrays for this filter and config."""
        self.filter = filter # filter for transmission
        # Default to full LW range if filter not recognized
        self.WRANGE = self.__WRANGE__.get(filter, np.array([2.4, 5.1]))
        self.disp_filter = self.__DISP_FILTER__.get(filter) # filter for dispersion
        self.config = config
        self.finalscale_drizzle = getattr(config, 'finalscale_drizzle', 0.035)
        self.pixfrac_drizzle    = getattr(config, 'pixfrac_drizzle',    0.8)

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
        """
        Load spectral-trace polynomial coefficients for ``fit_disp_order23``.

        The file encodes ``dy(x0, y0, dx)`` where dx = x_s − x0, dy = y_s − y0.

        Parameters
        ----------
        module:
            NIRCam module (``'A'`` or ``'B'``).
        grism:
            Grism identifier (``'R'`` or ``'C'``).

        Returns
        -------
        coeff_bestfit:
            1-D array of 30 best-fit trace polynomial coefficients.
        coeff_error:
            1-D array of 30 coefficient uncertainties.
        """
        path = os.path.join(
            self.config.cali_support_dir,
            "DISP_%s_mod%s_grism%s.dat" % (self.disp_filter, module, grism),
        )
        tb = ascii.read(path)
        return tb["col0"].data, tb["col1"].data
    
    def _load_displ(self, module, grism):
        """
        Load dispersion wavelength-solution coefficients for ``fit_disp_order32``.

        The file encodes ``dx(x0, y0, lambda_s)``.

        Parameters
        ----------
        module:
            NIRCam module (``'A'`` or ``'B'``).
        grism:
            Grism identifier (``'R'`` or ``'C'``).

        Returns
        -------
        coeff_bestfit:
            1-D array of 16 best-fit dispersion polynomial coefficients.
        coeff_error:
            1-D array of 16 coefficient uncertainties.
        """
        path = os.path.join(
            self.config.cali_support_dir,
            "DISPL_mod%s_grism%s.dat" % (module, grism),
        )
        tb = ascii.read(path)
        return tb["col0"].data, tb["col1"].data
    
    def _load_sens(self, module, grism):
        """
        Load the flux-sensitivity curve for a module/grism combination.

        Parameters
        ----------
        module:
            NIRCam module (``'A'`` or ``'B'``).
        grism:
            Grism identifier (``'R'`` or ``'C'``).

        Returns
        -------
        interpolate_sens:
            Callable that maps wavelength (µm) → sensitivity (DN/s/Jy).
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
        """
        Return the 16 dispersion polynomial coefficients for ``fit_disp_order32``.

        Parameters
        ----------
        module:
            NIRCam module (``'A'`` or ``'B'``).
        grism:
            Grism identifier (``'R'`` or ``'C'``).
        """
        idx = self.list_mod_pupil.index(f'{module}{grism}')
        return self.list_disp_coeff[idx]

    def get_trace_coeff(self, module, grism):
        """
        Return the 30 trace polynomial coefficients for ``fit_disp_order23``.

        Parameters
        ----------
        module:
            NIRCam module (``'A'`` or ``'B'``).
        grism:
            Grism identifier (``'R'`` or ``'C'``).
        """
        idx = self.list_mod_pupil.index(f'{module}{grism}')
        return self.list_trace_coeff[idx]
    
    def get_sensitivity(self, module, grism):
        """
        Return the sensitivity callable for the given module/grism combination.

        Parameters
        ----------
        module:
            NIRCam module (``'A'`` or ``'B'``).
        grism:
            Grism identifier (``'R'`` or ``'C'``).

        Returns
        -------
        interpolate_sens:
            Callable that maps wavelength (µm) → sensitivity (DN/s/Jy).
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
