"""
extraction.py
-------------
2-D spectral extraction, 1-D spectral extraction, and 2-D emission-line
cutout co-addition for NIRCam WFSS grism data.

Module layout (in order)
~~~~~~~~~~~~~~~~~~~~~~~~
Section 1 – Single-frame 2-D extraction
    ``extract_2d_spec``            – extract a 2-D spectrum from one grism frame.

Section 2 – Multi-frame 2-D spectral co-addition
    ``store_all_2d_spec``          – combine per-frame 2-D spectra into one FITS.
    ``resample_spec2d_wmin_wmax``  – resample a 2-D spectrum column into a wavelength bin.
    ``extract_2d_spec_worker``     – pool-safe driver: calls extract_2d_spec per source.

Section 3 – Direct-image cutout helpers and 1-D extraction
    ``_get_mosaic_cutout``         – cut a stamp from a local mosaic FITS file.
    ``_query_jades_cutout``        – stream a JADES HLSP stamp from MAST.
    ``extract_1d_spec_worker``     – pool-safe driver: 1-D extraction + diagnostic plot.

Section 4 – 2-D emission-line cutout extraction
    ``store_all_2d_emline``        – combine per-frame emission-line cutouts into FITS.
    ``_compute_grism_psf_frame``   – compute a stpsf PSF model for one grism frame.
    ``extract_2d_emline_worker_drizzle``– pool-safe driver: drizzle-coadd emission-line cutouts.
    ``extract_2d_emline_worker_simple`` - pool-safe driver: simple-coadd emission line cutouts without re-sampling. 

Comment about `GS_V3_PA`:
The `GS_V3_PA` equals the roll angle (northward from west) of the Sci-frame image for
non-obvious reasons. The definition of `GS_V3_PA` is the V3 axis position angle defined
eastward of north, measured at guide star region in NIRCam. Note that the Y axes of JWST
detectors generally are misaligned with the V3 axis, but NIRCam happens to sit at the 
center of JWST focal plane and its Y axis is only offset from V3 by ~0.02 degrees.
Therefore the roll angle (angle of +X from west to north) determination of this pipeline 
can not be applied to other JWST instruments without modification.
"""

from __future__ import annotations

import os
import time

import numpy as np

from astropy.io import fits, ascii
from astropy.table import Table, vstack, join
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.stats import SigmaClip, sigma_clipped_stats
from astropy.wcs import WCS
from astropy.modeling import models
from astropy.visualization import ZScaleInterval
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from scipy import ndimage, optimize, interpolate

from nircam_wfss.dispersion import grism_conf_preparation, get_position_ang_dispang_at_wave
from nircam_wfss.plotting import corner_text
from nircam_wfss.dispersion import GrismConf
from nircam_wfss.config import EML_LAB, NIRCAM_LW_PIXSCALE
from drizzle.resample import Drizzle

VALID_X_LEFT = 5
VALID_X_RIGHT = 2047 - 6
VALID_Y_BOTTOM = 5
VALID_Y_TOP = 2047 - 6


'''gaussian function'''
gauss = lambda x, x0, flux, FWHM : (flux / FWHM / 1.064467) * np.exp(-(x - x0)**2 / (2 * (FWHM/2.354820)**2))

# ---------------------------------------------------------------------------
# Single-frame 2-D extraction
# ---------------------------------------------------------------------------

def extract_2d_spec(
    img: np.ndarray,
    wrange: tuple[float, float],
    x0: float,
    y0: float,
    dxs: np.ndarray,
    dys: np.ndarray,
    wave: np.ndarray,
    img_wht: np.ndarray,
    img_dq: np.ndarray,
    header: fits.Header,
    pupil: str,
    module: str, 
    img_line: np.ndarray | None = None,
    aper: float = 10.0,
    unified_chirality: bool = True,
) -> fits.HDUList:
    """
    Extract a 2-D spectrum from one dispersed slitless spectroscopic image.

    The extraction uses a simple rectangular aperture centred on the
    wavelength trace.  Sub-pixel shifts along the cross-dispersion axis are
    corrected with a linear (order-1) ``scipy.ndimage.shift``.

    Comment on the output spectrum orientation or chirality: given the sophisticated
    design of JWST NIRCam WFSS and the methodology of this function, the output 
    spectra have a chirality. Defining the Z axis as the cross-product of the SCI 
    frame +Z = (+X) x (+Y), and putting the wavelength direction +W and the cross-
    dispersion direction +C on the X-Y plane, then the chirality is +1 if 
            [ (+W) x (+C) ] * (+Z) > 0,
    and -1 if otherwise. This doesn't matter for redshift measurement, but matters
    for morphology measurement, and kinematic lensing modeling in particular, which
    is sensitive to the direction of the dispersion, and requires consistent position 
    angle definition between grism and imaging, and also a convention consistent with 
    the WCS.

    Parameters
    ----------
    img:
        2-D grism science image.
    wrange:
        ``(wave_min, wave_max)`` in µm defining the extraction wavelength
        range.
    x0, y0:
        Source position in the direct image (pixels). Should be
        the position after astrometric correction and POM modulation.
    dxs, dys:
        Pixel offsets of the wavelength trace from (x0, y0), computed by
        :func:`~nircam_wfss.dispersion.grism_conf_preparation`.
    wave:
        Wavelength array (µm) associated with each (dx, dy) step.
    img_wht:
        Weight image (1/rms²).
    img_dq:
        Data-quality image (0 = good).
    header:
        Primary header of the grism FITS file (used to copy exposure
        metadata into the output).
    pupil: str
        Grism pupil: ``'R'`` (GRISMR, disperses along X) or ``'C'``
        (GRISMC, disperses along Y).
    module: str
        Detector module: ``'A'`` or ``'B'``.
    img_line: np.ndarray | None
        Continuum-subtracted emission-line image.  If ``None``, the LINE2D
        extension is filled with zeros.
    aper: float, optional
        Half-aperture width in the cross-dispersion direction (pixels).
        Defaults to 10, which corresponds to a 21-pixel wide aperture.
    unified_chirality: bool, optional
        If ``True``, the output spectra are flipped along the cross-dispersion 
        axis if the module-pupil combination is not AR, such that the output
        spectra is equivalent to AR up to a roll angle rotation of (0, 360) deg.
        If ``False``, the output spectra are not flipped such that the cross-
        dispersion direction is always +X or +Y in the original SCI frame. 
        For the ease of kinematic lensing modeling, we recommend set to True,
        which is the default behavior.
    Returns
    -------
    hdul:
        FITS ``HDUList`` with extensions:
        ``PRIMARY``, ``SPEC2D``, ``WHT2D``, ``DQ2D``, ``LINE2D``, ``WAVE``.

    Raises
    ------
    ValueError
        If fewer than 20 wavelength pixels fall within the requested range
        and detector boundaries.
    KeyError
        If ``pupil`` is not ``'R'`` or ``'C'``.
        If ``module`` is not ``'A'`` or ``'B'``.
    """
    if pupil not in ("C", "R"):
        raise KeyError('pupil must be "R" or "C"')
    if module not in ("A", "B"):
        raise KeyError('module must be "A" or "B"')

    w_min, w_max = wrange
    x_on_g = dxs + x0   # trace x coordinates on the grism image
    y_on_g = dys + y0   # trace y coordinates on the grism image
    aper_int = int(aper)

    # Valid pixel range (stay 5–6 px from the detector edge)
    if pupil == "R":
        args_eff = np.where(
            (wave >= w_min) & (wave <= w_max)
            & (x_on_g >= VALID_X_LEFT) & (x_on_g <= VALID_X_RIGHT)
            & (y_on_g >= VALID_Y_BOTTOM + aper_int) & (y_on_g <= VALID_Y_TOP - aper_int)
        )
    else:  # 'C'
        args_eff = np.where(
            (wave >= w_min) & (wave <= w_max)
            & (x_on_g >= VALID_X_LEFT + aper_int) & (x_on_g <= VALID_X_RIGHT - aper_int)
            & (y_on_g >= VALID_Y_BOTTOM) & (y_on_g <= VALID_Y_TOP)
        )

    if np.size(args_eff) <= 20:
        raise ValueError("Too few wavelength pixels in the valid detector region")

    n_wave = len(args_eff[0])
    n_cross = aper_int * 2 + 1

    tmp_spec_2d = np.zeros((n_wave, n_cross))
    tmp_wht_2d  = np.zeros((n_wave, n_cross))
    tmp_dq_2d   = np.zeros((n_wave, n_cross))
    tmp_line_2d = np.zeros((n_wave, n_cross))

    if img_line is None:
        img_line = np.zeros_like(img)

    for i, j in enumerate(args_eff[0]):
        # i: index of wavelength pixel in output arrays (tmp_spec_2d, etc.)
        # j: index of wavelength pixel in input arrays (dxs, dys, wave)
        if pupil == "R":
            # for each y-aperture slicing through dx trace, shift the y-aper along y 
            # to correct for sub-pixel trace position, then copy the shifted slice 
            # to output arrays. 
            tmp_x  = int(x_on_g[j])
            tmp_y1 = int(y_on_g[j] - aper_int - 1)
            tmp_y2 = int(y_on_g[j] + aper_int + 2)
            shift_sub = -(y_on_g[j] % 1)

            for arr in (img, img_wht, img_dq, img_line):
                arr.T[tmp_x, tmp_y1:tmp_y2] = ndimage.shift(
                    arr.T[tmp_x, tmp_y1:tmp_y2], shift_sub, order=1, mode="wrap"
                )
            # tmp_*_2d are of dimension (n_wave, n_cross), where n_cross direction is +Y
            tmp_spec_2d[i] = img.T[tmp_x][tmp_y1 + 1 : tmp_y2 - 1]
            tmp_wht_2d[i]  = img_wht.T[tmp_x][tmp_y1 + 1 : tmp_y2 - 1]
            tmp_dq_2d[i]   = img_dq.T[tmp_x][tmp_y1 + 1 : tmp_y2 - 1]
            tmp_line_2d[i] = img_line.T[tmp_x][tmp_y1 + 1 : tmp_y2 - 1]

        else:  # 'C'
            # for each x-aperture slicing through dy trace, shift the x-aper along x 
            # to correct for sub-pixel trace position, then copy the shifted slice 
            # to output arrays. 
            tmp_y  = int(y_on_g[j])
            tmp_x1 = int(x_on_g[j] - aper_int - 1)
            tmp_x2 = int(x_on_g[j] + aper_int + 2)
            shift_sub = -(x_on_g[j] % 1)

            for arr in (img, img_wht, img_dq, img_line):
                arr[tmp_y, tmp_x1:tmp_x2] = ndimage.shift(
                    arr[tmp_y, tmp_x1:tmp_x2], shift_sub, order=1, mode="wrap"
                )
            # tmp_*_2d are of dimension (n_wave, n_cross), where n_cross direction is +X
            tmp_spec_2d[i] = img[tmp_y, tmp_x1 + 1 : tmp_x2 - 1]
            tmp_wht_2d[i]  = img_wht[tmp_y, tmp_x1 + 1 : tmp_x2 - 1]
            tmp_dq_2d[i]   = img_dq[tmp_y, tmp_x1 + 1 : tmp_x2 - 1]
            tmp_line_2d[i] = img_line[tmp_y, tmp_x1 + 1 : tmp_x2 - 1]

    # Transpose to (cross-dispersion, wavelength) orientation
    # Direction of each dimensions of (cross-dispersion, wavelength) wrt SCI frame:
    # 1. AR: (cross-dispersion, wavelength) = (+Y, +X) (chirality +1)
    # 2. BR: (cross-dispersion, wavelength) = (+Y, -X) (chirality -1)
    # 3. AC: (cross-dispersion, wavelength) = (+X, +Y) (chirality -1)
    # 4. BC: (cross-dispersion, wavelength) = (+X, +Y) (chirality -1)
    tmp_spec_2d = tmp_spec_2d.T
    tmp_wht_2d  = tmp_wht_2d.T
    tmp_dq_2d   = tmp_dq_2d.T
    tmp_line_2d = tmp_line_2d.T
    chirality = 1 if (module, pupil) == ("A", "R") else -1
    if unified_chirality and (f"{module}{pupil}" != "AR"):
        # Now, we will flip the spectra along the cross-dispersion direction if
        # the chirality is -1, so that the output spectra chirality is equivalent
        # to GRISM AR up to a roll angle difference (0, 360) degrees.
        tmp_spec_2d = tmp_spec_2d[::-1]
        tmp_wht_2d  = tmp_wht_2d[::-1]
        tmp_dq_2d   = tmp_dq_2d[::-1]
        tmp_line_2d = tmp_line_2d[::-1]
        chirality = 1
        # And now the directions are
        # 1. AR: (cross-dispersion, wavelength) = (+Y, +X) (chirality +1)
        # 2. BR: (cross-dispersion, wavelength) = (-Y, -X) (chirality +1)
        # 3. AC: (cross-dispersion, wavelength) = (-X, +Y) (chirality +1)
        # 4. BC: (cross-dispersion, wavelength) = (-X, +Y) (chirality +1)

    # --- Build output FITS ---
    hdu = fits.PrimaryHDU()
    hdu.header["X0"]     = (np.float32(x0), "Reference position X in direct image")
    hdu.header["Y0"]     = (np.float32(y0), "Reference position Y in direct image")
    hdu.header["AUTHOR"] = ("Jiachuan Xu", "Author of this file")
    hdu.header["CHIRAL"] = (chirality, "Trace chirality: +1 or -1")
    hdu.header["TIME"]   = (
        time.strftime("%Y/%m/%d %H:%M:%S", time.localtime()), "Time of Creation"
    )
    meta_keys = [
        "TITLE", "PI_NAME", "CATEGORY", "SCICAT",
        "DATE-OBS", "TIME-OBS", "DATE-BEG", "DATE-END",
        "OBS_ID", "VISIT_ID", "PROGRAM", "OBSERVTN", "OBSLABEL",
        "OBSFOLDR", "GS_V3_PA",
        "EXPSTART", "EXPMID", "EXPEND",
        "READPATT", "NINTS", "NGROUPS", "EFFINTTM", "EFFEXPTM", "DURATION",
    ]
    for key in meta_keys:
        if key in header:
            hdu.header[key] = header.cards[key][1:]

    hdu_sci = fits.ImageHDU(np.float32(tmp_spec_2d), name="SPEC2D")
    hdu_sci.header["WAVE_1"]   = (float(wave[args_eff[0][0]]),
                                  "Wavelength (um) of first pixel")
    hdu_sci.header["D_WAVE"]   = (float(np.mean(np.diff(wave[args_eff[0]]))),
                                  "Wavelength step (um) per pixel")
    hdu_sci.header["COMMENTS"] = "wave = WAVE_1 + arange(0,NAXIS1) * D_WAVE"
    hdu_sci.header["APERTURE"] = (aper, "Aperture half-width in cross-dispersion (pix)")
    hdu_sci.header["PUPIL"]    = (pupil, "Pupil (R=GRISMR, C=GRISMC)")
    hdu_sci.header["MODULE"]   = (header.get("MODULE", "?"), "Detector module (A or B)")
    hdu_sci.header["DIFF_Y"]   = (float(y_on_g[args_eff[0]][0]), "Y_(full)-Y_(trim)")
    hdu_sci.header["DIFF_X"]   = (float(x_on_g[args_eff[0]][0]), "X_(full)-X_(trim)")

    hdu_wht = fits.ImageHDU(np.float32(tmp_wht_2d), name="WHT2D")
    hdu_dq  = fits.ImageHDU(np.int32(tmp_dq_2d),   name="DQ2D")

    hdu_line = fits.ImageHDU(np.float32(tmp_line_2d), name="LINE2D")
    hdu_line.header["COMMENTS"] = "Extracted from continuum-subtracted map"

    tb_wave = Table(
        data=[
            wave[args_eff],
            x_on_g[args_eff],
            y_on_g[args_eff],
            dxs[args_eff],
            dys[args_eff],
        ],
        names=["wavelength", "xs", "ys", "dxs", "dys"],
    )
    tb_wave["wavelength"].info.format = ".6f"
    for col in tb_wave.colnames[1:]:
        tb_wave[col].info.format = ".3f"

    hdul = fits.HDUList([hdu, hdu_sci, hdu_wht, hdu_dq, hdu_line,
                         fits.BinTableHDU(tb_wave, name="WAVE")])
    return hdul


# ---------------------------------------------------------------------------
# Multi-frame co-addition into a single output FITS
# ---------------------------------------------------------------------------

def store_all_2d_spec(
    fits_list: list[fits.HDUList],
    pupils: list[str],
    modules: list[str],
    paths: list[str],
    output: str = "spec_2d.fits",
    coord=None,
    grism_filter: str | None = None,
    info_table=None,
    overwrite: bool = False,
) -> fits.HDUList:
    """
    Combine per-exposure 2-D spectral extractions into one FITS file.
    Each input spectrum FITS file contains the extensions: ``PRIMARY``, 
    ``SPEC2D``, ``WHT2D``, ``DQ2D``, ``LINE2D`` (optional), and ``WAVE``.
    All the extensions except the primary header are concatenated into 
    the output FITS with new extension names like ``SPEC2D-0``, ``SPEC2D-1``, 
    etc. The output file is written to `output` when ``overwrite=True``. 

    Parameters
    ----------
    fits_list:
        List of ``HDUList`` objects returned by :func:`extract_2d_spec`.
    pupils:
        Pupil label (``'R'`` or ``'C'``) for each entry in ``fits_list``.
    modules:
        Module label (``'A'`` or ``'B'``) for each entry.
    paths:
        Path to the source grism FITS file for each entry.
    output:
        Output file path.
    coord:
        ``astropy.coordinates.SkyCoord`` of the source (optional).
    grism_filter:
        Filter name (e.g. ``'F444W'``) to record in the header.
    info_table:
        A single-row ``astropy.table.Row`` whose columns are written to the
        primary header under ``HIERARCH`` keywords.
    overwrite:
        Overwrite the output file if it exists. Default is False.

    Returns
    -------
    ind_hdul:
        The assembled ``HDUList`` (also written to disk when ``overwrite=True``).
    """
    ind_hdul = fits.HDUList([fits.PrimaryHDU()])

    for l, x in enumerate(fits_list):
        # Primary header / SPEC2D / WHT2D / DQ2D / [LINE2D] / WAVE
        ind_hdul.append(x[1])
        ind_hdul[-1].header["EXTNAME"]  = "SPEC2D-%d" % l
        for card in ("x0", "y0", "CHIRAL"):
            ind_hdul[-1].header[card] = x[0].header[card]
        ind_hdul[-1].header["PUPIL"]    = pupils[l]
        ind_hdul[-1].header["MODULE"]   = modules[l]
        ind_hdul[-1].header["DATAPATH"] = os.path.basename(paths[l])

        ind_hdul.append(x[2])
        ind_hdul[-1].header["EXTNAME"] = "WHT2D-%d" % l
        ind_hdul.append(x[3])
        ind_hdul[-1].header["EXTNAME"] = "DQ2D-%d" % l

        # LINE2D is extension 4 when present; WAVE is the last extension
        if len(x) == 5:
            # No LINE2D – last extension is WAVE
            ind_hdul.append(x[4])
            ind_hdul[-1].header["EXTNAME"] = "WAVE-%d" % l
        else:
            ind_hdul.append(x[4])
            ind_hdul[-1].header["EXTNAME"] = "LINE2D-%d" % l
            ind_hdul.append(x[5])
            ind_hdul[-1].header["EXTNAME"] = "WAVE-%d" % l

    # Summary table
    stats_table = Table(
        names=["id", "pupil", "module", "datapath"],
        data=[
            list(range(len(fits_list))),
            pupils,
            modules,
            [os.path.basename(p) for p in paths],
        ],
    )
    ind_hdul.append(fits.BinTableHDU(stats_table, name="STATS"))

    # Primary header metadata
    ind_hdul[0].header["DIRNAME"] = (
        os.path.dirname(paths[0]), "Directory of original grism data"
    )
    if coord is not None:
        ind_hdul[0].header["RA0"]  = (float(coord.ra.deg),  "Source RA (deg)")
        ind_hdul[0].header["DEC0"] = (float(coord.dec.deg), "Source Dec (deg)")
    ind_hdul[0].header["N_COADD"] = (len(pupils), "Total coadded frames")
    ind_hdul[0].header["N_R"] = (sum(p == "R" for p in pupils), "Frames from GRISMR")
    ind_hdul[0].header["N_C"] = (sum(p == "C" for p in pupils), "Frames from GRISMC")
    ind_hdul[0].header["AUTHOR"] = ("Jiachuan Xu", "Author")
    ind_hdul[0].header["TIME"]   = (
        time.strftime("%Y/%m/%d %H:%M:%S", time.localtime()), "Creation time"
    )
    if grism_filter is not None:
        ind_hdul[0].header["FILTER"] = (grism_filter, "Filter name")

    if info_table is not None:
        ind_hdul[0].header["COMMENTS"] = "Catalog information below:"
        for col in info_table.colnames:
            val = info_table[col]
            if isinstance(val, np.ma.core.MaskedConstant):
                continue
            try:
                if isinstance(val, (str, np.str_)):
                    ind_hdul[0].header["HIERARCH " + col] = str(val)
                elif np.isnan(val):
                    ind_hdul[0].header["HIERARCH " + col] = "nan"
                elif np.isinf(val):
                    ind_hdul[0].header["HIERARCH " + col] = "inf"
                else:
                    ind_hdul[0].header["HIERARCH " + col] = val
            except Exception:
                pass

    if overwrite:
        ind_hdul.writeto(output, overwrite=True)

    return ind_hdul


# ---------------------------------------------------------------------------
# Wavelength-bin resampling (for co-addition at uniform wavelength grid)
# ---------------------------------------------------------------------------

def resample_spec2d_wmin_wmax(
    x: int,
    i: int,
    tmp_wave_2d: list[np.ndarray],
    tmp_ind_fits_list: list[fits.HDUList],
    wave_sample: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Resample a single wavelength column of a 2-D spectrum.

    This function is designed to be called via ``Pool.starmap`` or in a
    simple loop.  It replaces the version in the original notebook that
    used global variables.

    .. note::
        This function is not called anywhere in the current codebase.
        ``extract_2d_spec_worker`` inlines equivalent logic directly.

    Parameters
    ----------
    x:
        Index into ``tmp_ind_fits_list`` (which grism exposure to use).
    i:
        Index into ``wave_sample`` (which wavelength bin to fill).
    tmp_wave_2d:
        List of 1-D wavelength arrays, one per exposure.
    tmp_ind_fits_list:
        List of ``HDUList`` objects, each from :func:`store_all_2d_spec`.
    wave_sample:
        Uniformly-spaced output wavelength grid (edges).

    Returns
    -------
    tmp_spec_w_1:
        Weighted-mean spatial slice of SCI for this wavelength bin.
    tmp_wht_w_1:
        Summed weight for this bin.
    tmp_cov_w_1:
        Integer coverage (1 where covered, 0 elsewhere).
    tmp_line_w_1:
        Weighted-mean spatial slice of LINE2D for this bin.
    """
    tmp_w_min = wave_sample[i]
    tmp_w_max = wave_sample[i + 1]

    in_bin = (tmp_wave_2d[x] > tmp_w_min) & (tmp_wave_2d[x] <= tmp_w_max)
    arg_in = tuple([in_bin])

    tmp_spec_w = tmp_ind_fits_list[x][1].data.T[arg_in]
    tmp_wht_w  = np.nan_to_num(
        tmp_ind_fits_list[x][2].data.T[arg_in], posinf=0, neginf=0
    )
    tmp_dq_w   = tmp_ind_fits_list[x][3].data.T[arg_in]

    # Loose DQ flag: only mask "do not use" (bit 0)
    tmp_wht_w[tmp_dq_w % 2 == 1] = 0

    # LINE2D: extension 4 when no LINE2D, else 5
    if len(tmp_ind_fits_list[x]) == 5:
        tmp_line_w = tmp_ind_fits_list[x][1].data.T[arg_in]
    else:
        tmp_line_w = tmp_ind_fits_list[x][4].data.T[arg_in]

    wht_sum = np.nansum(tmp_wht_w, axis=0)
    tmp_spec_w_1 = np.nansum(tmp_spec_w * tmp_wht_w, axis=0) / wht_sum
    tmp_wht_w_1  = wht_sum
    tmp_cov_w_1  = np.int8(wht_sum != 0)
    tmp_line_w_1 = np.nansum(tmp_line_w * tmp_wht_w, axis=0) / wht_sum

    return tmp_spec_w_1, tmp_wht_w_1, tmp_cov_w_1, tmp_line_w_1



def extract_2d_spec_worker(
    grism_idx_per_source: list[str],
    POM_path_per_source: list[str],
    all_v1p5_list: list[str],
    source_item: Table.Row,
    grism_conf: "GrismConf",
    aper: float,
    filter: str,
    extraction_dir: str,
    bunit: str,
    unified_chirality: bool = True,
    overwrite: bool = False,
) -> None:
    """
    Extract 2-D spectra for a single source.

    This function calls :func:`extract_2d_spec` to get the 2-D spectrum, weight map, 
    DQ map, and optionally the emission-line map for each grism frame where the source 
    is observable, and then saves each spectrum as a FITS file in `extraction_dir`,
    both as a per-frame compilation and co-adds grouped by module-pupil combination. 
    This function is designed to be called via ``Pool.starmap`` or in a simple loop.

    Parameters
    ----------
    grism_idx_per_source:
        List of grism frame indices where the source is observable. The index can be
        used to find the corresponding grism frame in `all_v1p5_list`.
    POM_path_per_source:
        List of paths to the POM catalogs where the source is observable. The length
        of this list should be the same as `grism_idx_per_source`.
    all_v1p5_list:
        List of paths to all v1.5 grism FITS files (used to find the corresponding 
        grism frame for each POM catalog).
    source_item:
        A single row from the source catalog table, containing at least
        'ID', 'RA', and 'DEC' columns.
    grism_conf: `GrismConf`
        Grism configuration object containing dispersion information.
    aper: float
        Aperture size for the extraction.
    filter: str
        Filter name for the extraction.
    extraction_dir: str
        Directory where the extracted spectra will be saved.
    bunit: str
        Brightness unit of the extracted spectra ('mJy' or 'DN/s').
    unified_chirality: bool
        Whether to unify the chirality of the output spectra across different 
        module-pupil combinations. See the `unified_chirality` parameter in 
        `extract_2d_spec` for details. Default is True.
    overwrite: bool
        Whether to overwrite existing 2D spectra files. Default is False.
    """
    if len(grism_idx_per_source) == 0: 
        print(' >> no spec found; skipping source %s' % source_item['ID'])
        return
    assert len(grism_idx_per_source) == len(POM_path_per_source), \
        "Length of grism_idx_per_source and POM_path_per_source should be the same."
    
    source_coord = SkyCoord(source_item["RA"], source_item["DEC"], unit=(u.deg, u.deg))
    source_id = source_item["ID"]


    spec2d_list     = [] # extracted 2D spectrum (HDUList) per frame
    fits_name_list  = [] # v1.5 filename per frame
    module_list     = [] # module per frame (A or B)
    pupil_list      = [] # pupil per frame (R or C)

    ''' 2D spectrum extraction & compilation
    In this step, we will loop through each grism exposure frames and extract 2D spectrum 
    for the source if the source is observable in that frame. Depending on the unit required,
    the extracted 2D spectra may or may not be corrected for sensitivity. Each temporary 2D 
    spectrum from each exposure will be saved in `spec2d_list`, and the corresponding grism 
    filename, module, and pupil information will be saved in `fits_name_list`, `module_list`, 
    and `pupil_list` respectively.
    Note that each item in `spec2d_list` is a HDUList containing the 2D spectrum (SPEC2D), 
    weight map (WHT2D), DQ map (DQ2D), and optionally the continuum-subtracted emission-line 
    map (LINE2D) and wavelength table (WAVE) for that frame.
    '''
    # For each grism frame where the source is observable
    for j, POM_fn in enumerate(POM_path_per_source):
        # get the reduced v1.5 grism filename
        grism_fn = all_v1p5_list[grism_idx_per_source[j]]
        # get the POM catalog for detector position in direct imaging
        POM_cat = ascii.read(POM_fn)

        ### Spectrum and Direct Image Position Retrieval
        # Read grism data (w/ and w/o continuum subtraction) & header
        image = fits.getdata(grism_fn, "sci") # grism w/ continuum
        try:
            emline = fits.getdata(grism_fn, 'emline') # Line-only Grism SCI image data
        except KeyError:
            emline = image
        data_quality = fits.getdata(grism_fn, "dq")
        primary_hd = fits.getheader(grism_fn)
        sci_hd = fits.getheader(grism_fn, "sci")
        _filter, module, pupil = primary_hd["FILTER"], primary_hd["MODULE"], primary_hd["PUPIL"][-1]
        if _filter != filter:
            print(' >> filter mismatch; skipping frame %s' % grism_fn)
            continue
        # Read or Generate Grism WHT data:
        weight_path = grism_fn.replace('lv1.5.fits', 'wht.fits')
        if os.path.isfile(weight_path):
            weight = fits.getdata(weight_path)
        else:
            ### 2D grism WHT file not found, generate from Error map:
            weight = fits.getdata(grism_fn, 'err')
            weight[weight == 0] = np.nan
            weight = weight**-2
            weight_hdul = fits.HDUList([fits.PrimaryHDU(), fits.ImageHDU(weight, name = 'WHT')])
            weight_hdul[1].header = sci_hd.copy()
            # Write atomically via temp file to avoid multiprocessing race conditions:
            # multiple workers share the same grism file and may all try to write
            # the weight file simultaneously; os.rename() is atomic on POSIX filesystems.
            tmp_path = weight_path + '.tmp.%d' % os.getpid()
            weight_hdul.writeto(tmp_path, overwrite=True)
            weight_hdul.close()
            try:
                os.rename(tmp_path, weight_path)
            except Exception:
                os.remove(tmp_path)  # another worker finished first; clean up our temp
            weight = fits.getdata(weight_path)
        # direct image position 
        item_POM = POM_cat[POM_cat['Index'] == source_id][0]
        x0 = item_POM['pixel_x']
        y0 = item_POM['pixel_y']
        print('%s(%.1f, %.1f) ' % (pupil, x0, y0), end = ' ')

        ### Wavelength solution
        # spectral tracing parameters
        disp_coeff = grism_conf.get_disp_coeff(module, pupil)
        trace_coeff = grism_conf.get_trace_coeff(module, pupil)
        dxs, dys, wavs = grism_conf_preparation(x0 = x0,  y0 = y0, pupil = pupil, 
                        fit_opt_fit = trace_coeff, w_opt = disp_coeff)
        # correct wavelength to barycentric coordinate
        wavs = (1 + sci_hd['velosys']/299792458) * wavs

        ### Extract 2D spectrum on each single frame
        try:
            tmp_spec_2D = extract_2d_spec(
                img = image, wrange = grism_conf.WRANGE, 
                x0 = x0, y0 = y0, dxs = dxs, dys = dys, wave = wavs, 
                img_wht = weight, img_dq = data_quality,
                img_line = emline, ## add EMLINE extension
                header = primary_hd, pupil = pupil, module = module, 
                aper = aper, unified_chirality = unified_chirality)
        except ValueError:
            continue
        spec2d_list.append(tmp_spec_2D)
        module_list.append(module)
        fits_name_list.append(grism_fn)
        pupil_list.append(pupil)

    ''' Save 2D spectra of all individual extraction 
    Note that this step only save the extracted 2D spectra from each frame in DN/s.
    '''
    if len(spec2d_list) == 0: 
        print(' >> no spec found; ')
        return
    else: 
        print(' >> N=%d, x0=%.1f, y0=%.1f' % (len(fits_name_list), x0, y0), end = '')
        specs2d_compile_name = os.path.join(extraction_dir, 
                                            'allspec_2d_%s_ID%s.fits' % (filter, source_id))
        store_all_2d_spec(fits_list = spec2d_list, pupils = pupil_list, 
                            modules = module_list, paths = fits_name_list, 
                            output = specs2d_compile_name, coord = source_coord, 
                            grism_filter = filter,  info_table = source_item, 
                            overwrite = overwrite)
        print(' >> save all extracted spec2d ', end = '')

    ''' Coadd 2D spectra from each frame by module-pupil separately 
    Note that depending on the unit required, the coadded 2D spectra may or may not be 
    corrected for sensitivity. If the unit is in mJy, the coadded 2D spectra will be 
    corrected for sensitivity; if the unit is in DN/s, the coadded 2D spectra will not 
    be corrected for sensitivity. The final coadded 2D spectra will be saved in 
    `spec_2d_{filter}_ID{source_id}_{module-pupil|all}coadd.fits`
    '''
    # Correct for sensitivity (depending on the unit required). 
    for k, tmp_fits in enumerate(spec2d_list[:]):
        # default sensitivity: DN/s per Jy
        if len(tmp_fits) == 5:
            tmp_wavelength = tmp_fits[4].data['wavelength']
        else: # if has a line map
            tmp_wavelength = tmp_fits[5].data['wavelength']
        tmp_f_sens = grism_conf.get_sensitivity(tmp_fits[1].header['MODULE'], 
                                                tmp_fits[1].header['PUPIL'])
        if bunit == 'mJy':
            f_sens = tmp_f_sens(tmp_wavelength) * 1.0e-3 # DN/s per mJy
        elif bunit == 'DN/s':
            f_sens = np.ones_like(tmp_wavelength)
        else:
            raise ValueError('bunit must be "mJy" or "DN/s"')
        # change units in spec2d extension
        spec2d_list[k][1].header['BUNIT'] = (bunit, 'Brightness Unit')
        spec2d_list[k][1].data = tmp_fits[1].data / f_sens
        spec2d_list[k][2].data = tmp_fits[2].data * f_sens**2
        if len(tmp_fits) != 5:
            # change units in emline extension
            spec2d_list[k][4].header['BUNIT'] = (bunit, 'Brightness Unit')
            spec2d_list[k][4].data = tmp_fits[4].data / f_sens
    # combine 2D spectra of all frames into one common grid
    d_wave = 0.001
    wave_sample = np.arange(grism_conf.WRANGE[0], grism_conf.WRANGE[1] + d_wave, d_wave)
    wave_sample_c = wave_sample[:-1] + d_wave / 2.
    N_common_wave = len(wave_sample) - 1
    N_spatial = len(spec2d_list[0][1].data)
    spec2d_list = np.array(spec2d_list, dtype=object)
    fits_name_list = np.array(fits_name_list)
    module_list = np.array(module_list)
    pupil_list = np.array(pupil_list)

    # stacked by module-pupil separately
    for module in ['A', 'B']:
        for pupil in ['R', 'C']:
            print(' >> coadd spec2d of %s%s ' % (module, pupil), end = '')
            # which exposures in `spec2d_list` belong to this module-pupil combination?
            ind_spec2d_list = np.array([i for i, x in enumerate(spec2d_list) \
                if x[1].header['MODULE'] == module and x[1].header['PUPIL'] == pupil], dtype=int)
            if len(ind_spec2d_list) == 0: 
                print(' >> no spec found; ')
                continue
            else: 
                N_exp = len(ind_spec2d_list)
                print(' >> N=%d, ' % N_exp, end = '')
                arr_spec_2d = np.zeros((N_common_wave, N_exp, N_spatial))
                arr_wht_2d  = np.zeros_like(arr_spec_2d)
                arr_cov_2d  = np.zeros_like(arr_spec_2d)
                arr_line_2d = np.zeros_like(arr_spec_2d)
                # for each wavelength bin in the common grid, deposit the native spectrum
                # by picking the native pixels that fall into this wavelength bin and taking 
                # their weighted average.
                # Note that this step does NOT re-sample!!!
                for local_iexp, iexp in enumerate(ind_spec2d_list):
                    native_spec = spec2d_list[iexp][1].data
                    native_wave = spec2d_list[iexp][-1].data['wavelength']
                    native_weight = spec2d_list[iexp][2].data
                    native_dq = spec2d_list[iexp][3].data
                    if len(spec2d_list[iexp]) == 5:
                        native_line = spec2d_list[iexp][1].data
                    else:
                        native_line = spec2d_list[iexp][4].data
                    for iwave in range(N_common_wave):
                        tmp_w_min = wave_sample[iwave]
                        tmp_w_max = wave_sample[iwave + 1]
                        # Does NOT re-sample!!!
                        in_bin = (native_wave > tmp_w_min) & (native_wave <= tmp_w_max)
                        arg_in = tuple([in_bin])
                        tmp_spec_w = native_spec.T[arg_in]
                        tmp_wht_w  = np.nan_to_num(
                            native_weight.T[arg_in], posinf=0, neginf=0
                        )
                        tmp_dq_w   = native_dq.T[arg_in]
                        # Loose DQ flag: only mask "do not use" (bit 0)
                        tmp_wht_w[tmp_dq_w % 2 == 1] = 0
                        tmp_line_w = native_line.T[arg_in]
                        wht_sum = np.nansum(tmp_wht_w, axis=0)
                        arr_spec_2d[iwave, local_iexp] = np.nansum(tmp_spec_w * tmp_wht_w, axis=0) / wht_sum
                        arr_wht_2d[iwave, local_iexp]  = wht_sum
                        arr_cov_2d[iwave, local_iexp]  = np.int8(wht_sum != 0)
                        arr_line_2d[iwave, local_iexp] = np.nansum(tmp_line_w * tmp_wht_w, axis=0) / wht_sum
                # the coadded spectra are sigma-clipped and weighted averaged across the 
                # exposures in each wavelength bin.
                sigma_clip = SigmaClip(sigma = 2.5, maxiters = 5)
                if arr_spec_2d.shape[1] > 2:
                    arr_spec_2d = sigma_clip(arr_spec_2d, axis = 1, masked = False)
                    arr_line_2d = sigma_clip(arr_line_2d, axis = 1, masked = False)
                arr_wht_2d[np.where(np.isnan(arr_spec_2d))] = 0
                arr_wht_2d[np.where(np.isnan(arr_line_2d))] = 0
                # weighted mean 2d spectra:
                coadd_spec_2d = np.nansum(arr_spec_2d * arr_wht_2d, axis = 1) / np.nansum(arr_wht_2d, axis = 1) # sum_arr_wht_2d # 
                coadd_wht_2d = np.nansum(arr_wht_2d, axis = 1)
                coadd_cov_2d = np.nansum(arr_cov_2d, axis = 1)
                coadd_line_2d = np.nansum(arr_line_2d * arr_wht_2d, axis = 1) / np.nansum(arr_wht_2d, axis = 1) # sum_arr_wht_2d # 
                # save coadded 2D spectrum of this module-pupil combination
                tmp_tb_cov = Table(
                    names = ['index', 'name', 'x0', 'y0', 'module', 'pupil', 'DIFF_X', 'DIFF_Y', 
                             'wave_0', 'wave_1', 'EFFEXPTM', 'GS_V3_PA', 'CHIRAL', 'DATE-BEG', 'DATE-END'],
                    data = [ind_spec2d_list,
                            [tmp_grism.split('/')[-1][:-16] for tmp_grism in fits_name_list[ind_spec2d_list]],
                            [x[0].header['x0'] for x in spec2d_list[ind_spec2d_list]],
                            [x[0].header['y0'] for x in spec2d_list[ind_spec2d_list]],
                            module_list[ind_spec2d_list],
                            pupil_list[ind_spec2d_list],
                            [x[1].header['DIFF_X'] for x in spec2d_list[ind_spec2d_list]],
                            [x[1].header['DIFF_Y'] for x in spec2d_list[ind_spec2d_list]],
                            [np.round(wave_sample_c[np.sum(arr_wht_2d, axis=-1)[:, n_] != 0][0], 4)
                                if np.any(np.sum(arr_wht_2d, axis=-1)[:, n_] != 0) else np.nan
                                for n_ in range(N_exp)],
                            [np.round(wave_sample_c[np.sum(arr_wht_2d, axis=-1)[:, n_] != 0][-1], 4)
                                if np.any(np.sum(arr_wht_2d, axis=-1)[:, n_] != 0) else np.nan
                                for n_ in range(N_exp)],
                            [x[0].header['EFFEXPTM'] for x in spec2d_list[ind_spec2d_list]],
                            [x[0].header['GS_V3_PA'] for x in spec2d_list[ind_spec2d_list]],
                            [x[0].header['CHIRAL'] for x in spec2d_list[ind_spec2d_list]],
                            [x[0].header['DATE-BEG'] for x in spec2d_list[ind_spec2d_list]],
                            [x[0].header['DATE-END'] for x in spec2d_list[ind_spec2d_list]],
                        ]
                )
                hdu = fits.PrimaryHDU()
                hdu.header['ID']     = (source_id,           'Source ID')
                hdu.header['RA0']    = (source_coord.ra.value,  'Reference position RA in direct image')
                hdu.header['DEC0']   = (source_coord.dec.value, 'Reference position Dec in direct image')
                hdu.header['N_coadd'] = (N_exp, 'Number of coadded frames')
                hdu.header['EFFEXPTM'] = (np.nansum(tmp_tb_cov['EFFEXPTM']), 'Maximum exposure time [s]')
                hdu.header['author'] = ("Jiachuan Xu", 'Author of this file')
                hdu.header['time'] = (time.strftime("%Y/%m/%d %H:%M:%S",  time.localtime()), 'Time of Creation')
                hdu.header['filter'] = (filter, 'Filter name')
                hdu.header['module'] = (module, 'Detector module (A or B)')
                hdu.header['pupil'] = (pupil, 'Pupil (R=GRISMR, C=GRISMC)')
                for cardname in ['PROGRAM', 'OBSERVTN', 'OBSLABEL', 'OBSFOLDR', 'GS_V3_PA']:
                    if cardname in primary_hd: hdu.header[cardname] = primary_hd.cards[cardname][1:]
                hdu.header['GS_V3_PA'] = (np.nanmean(tmp_tb_cov['GS_V3_PA']), 'V3 position angle of the grism pointing (deg)')
                
                ### copy source catalog information to primary header
                hdu.header['COMMENTS'] = 'Belows are information taken from input catalog:'
                for x in source_item.colnames:
                    if type(source_item[x]) == np.ma.core.MaskedConstant: continue
                    elif type(source_item[x]) != np.str_ :
                        try:
                            val = float(source_item[x])
                            if np.isnan(val) : hdu.header['HIERARCH ' + x] = 'nan'
                            elif np.isinf(val) : hdu.header['HIERARCH ' + x] = 'inf'
                            else: hdu.header['HIERARCH ' + x] = source_item[x]
                        except (TypeError, ValueError):
                            hdu.header['HIERARCH ' + x] = source_item[x]
                    else: hdu.header['HIERARCH ' + x] = source_item[x]
                ### save coadded data
                # image extension
                hdu_sci = fits.ImageHDU(coadd_spec_2d.T, name = 'SPEC2D')
                hdu_sci.header['wave_1'] = (wave_sample_c[0], 'Wavelength (um) of first pixel')
                hdu_sci.header['d_wave'] = (wave_sample_c[1] - wave_sample_c[0], 'Wavelength Difference (um) between each pixel')
                hdu_sci.header['comments'] = ('wave = wave_1 + np.arange(0, NAXIS1, 1) * d_wave')
                hdu_sci.header['pixscale'] = (NIRCAM_LW_PIXSCALE, 'Pixel scale in undispersed direction (arcsec)')
                hdu_sci.header['aperture'] = (aper, 'Aperture radius in undispersed direction (pixel)')
                hdu_sci.header['BUNIT'] = (bunit, 'Brightness Unit')
                # weight extension
                hdu_wht = fits.ImageHDU(coadd_wht_2d.T, name = 'WHT2D')
                hdu_wht.header['comments'] = ('Weight image; ERR = WHT^(-0.5)')
                # number of coverage per pixel extension
                hdu_cov = fits.ImageHDU(coadd_cov_2d.T, name = 'COV2D')
                hdu_cov.header['comments'] = ('Coverage image')
                # line-only SCI extension
                hdu_line = hdu_sci.copy()
                hdu_line.data = coadd_line_2d.T
                hdu_line.header['extname'] = 'LINE2D'
                hdu_line.header['comments'] = ('Line-only image extracted on continuum-filtered 2D data')
                hdu_line.header['BUNIT'] = (bunit, 'Brightness Unit')
                # statistics tab
                hdu_tab = fits.BinTableHDU(tmp_tb_cov, name = 'STATS')
                hdu_tab.header['COMMENT'] = 'name:     name of simulated image'
                hdu_tab.header['COMMENT'] = 'x0 / y0:  reference position (i.e., in direct image)'
                hdu_tab.header['COMMENT'] = '   (in reality this is registered to the wcs of grism image)'
                hdu_tab.header['COMMENT'] = 'DIFF_X:   X_(full)-X_(trim)'
                hdu_tab.header['COMMENT'] = 'DIFF_Y:   Y_(full)-Y_(trim)'
                hdu_tab.header['COMMENT'] = 'wave_0:   minimum wavelength (micron) in this coverage'
                hdu_tab.header['COMMENT'] = 'wave_1:   maximum wavelength (micron) in this coverage'
                hdul = fits.HDUList([hdu, hdu_sci, hdu_wht, hdu_cov, hdu_line, hdu_tab])
                # save into file
                sub_coadd_fn = os.path.join(extraction_dir,
                    'spec_2d_%s_ID%s_%scoadd.fits' % (filter, source_id, f'{module}{pupil}')
                )
                hdul.writeto(sub_coadd_fn, overwrite = True)
                print(' >> save module-pupil stacked spec2d; %s' % sub_coadd_fn)

    ''' Coadd 2D spectra from each frame all combined '''
    primary_hd = None
    all_spec_2d = None
    all_wht_2d = None
    all_line_2d = None
    all_table = None
    for module in ['A', 'B']:
        for pupil in ['R', 'C']:
            sub_coadd_fn = os.path.join(extraction_dir,
                    'spec_2d_%s_ID%s_%scoadd.fits' % (filter, source_id, f'{module}{pupil}')
                )
            if os.path.isfile(sub_coadd_fn):
                with fits.open(sub_coadd_fn) as tmp_hdul:
                    if primary_hd is None:
                        primary_hd = tmp_hdul[0].header.copy()
                    if all_spec_2d is None:
                        all_wht_2d = tmp_hdul[2].data
                        all_spec_2d = tmp_hdul[1].data * tmp_hdul[2].data
                        all_line_2d = tmp_hdul[4].data * tmp_hdul[2].data
                        all_table = Table(tmp_hdul[5].data)
                    else:
                        all_wht_2d = np.nansum([all_wht_2d, tmp_hdul[2].data], axis=0)
                        all_spec_2d = np.nansum([all_spec_2d, 
                                                 tmp_hdul[1].data * tmp_hdul[2].data], axis=0)
                        all_line_2d = np.nansum([all_line_2d, 
                                                 tmp_hdul[4].data * tmp_hdul[2].data], axis=0)
                        all_table = vstack([all_table, Table(tmp_hdul[5].data)])
    all_cov_2d = np.int8(all_wht_2d != 0)
    all_spec_2d = all_spec_2d / all_wht_2d
    all_line_2d = all_line_2d / all_wht_2d
    primary_hd["N_coadd"] = (len(all_table), 'Number of coadded frames')
    primary_hd["EFFEXPTM"] = (np.nansum(all_table['EFFEXPTM']), 'Maximum exposure time [s]')
    primary_hd["GS_V3_PA"] = (np.nanmean(all_table['GS_V3_PA']), 'V3 position angle of the grism pointing (deg)')
    primary_hd["N_R"] = (sum(all_table['pupil'] == 'R'), 'Number of frames from GRISMR')
    primary_hd["N_C"] = (sum(all_table['pupil'] == 'C'), 'Number of frames from GRISMC')
    hdu = fits.PrimaryHDU(header = primary_hd)
    # science extension
    hdu_sci = fits.ImageHDU(all_spec_2d, name = 'SPEC2D')
    hdu_sci.header['wave_1'] = (wave_sample_c[0], 'Wavelength (um) of first pixel')
    hdu_sci.header['d_wave'] = (wave_sample_c[1] - wave_sample_c[0], 'Wavelength Difference (um) between each pixel')
    hdu_sci.header['comments'] = ('wave = wave_1 + np.arange(0, NAXIS1, 1) * d_wave')
    hdu_sci.header["BUNIT"] = (bunit, 'Brightness Unit')
    # weight extension
    hdu_wht = fits.ImageHDU(all_wht_2d, name = 'WHT2D')
    hdu_wht.header['comments'] = ('Weight image; ERR = WHT^(-0.5)')
    # number of coverage per pixel extension
    hdu_cov = fits.ImageHDU(all_cov_2d, name = 'COV2D')
    hdu_cov.header['comments'] = ('Coverage image')
    # line-only SCI extension
    hdu_line = hdu_sci.copy()
    hdu_line.data = all_line_2d
    hdu_line.header['extname'] = 'LINE2D'
    hdu_line.header['comments'] = ('Line-only image extracted on continuum-filtered 2D data')
    hdu_line.header['BUNIT'] = (bunit, 'Brightness Unit')
    # statistics tab
    hdu_tab = fits.BinTableHDU(all_table, name = 'STATS')
    hdu_tab.header['COMMENT'] = 'name:     name of simulated image'
    hdu_tab.header['COMMENT'] = 'x0 / y0:  reference position (i.e., in direct image)'
    hdu_tab.header['COMMENT'] = '   (in reality this is registered to the wcs of grism image)'
    hdu_tab.header['COMMENT'] = 'DIFF_X:   X_(full)-X_(trim)'
    hdu_tab.header['COMMENT'] = 'DIFF_Y:   Y_(full)-Y_(trim)'
    hdu_tab.header['COMMENT'] = 'wave_0:   minimum wavelength (micron) in this coverage'
    hdu_tab.header['COMMENT'] = 'wave_1:   maximum wavelength (micron) in this coverage'
    hdu_tab.header['COMMENT'] = 'EFFEXPTM: effective exposure time (s)'
    hdu_tab.header['COMMENT'] = 'GS_V3_PA: V3 position angle of the grism pointing (deg)'
    hdul = fits.HDUList([hdu, hdu_sci, hdu_wht, hdu_cov, hdu_line, hdu_tab])
    # save all coadds into file
    coadd_name = os.path.join(extraction_dir, 'spec_2d_%s_ID%s_allcoadd.fits' % (filter, source_id))
    hdul.writeto(coadd_name, overwrite = True)
    print(' >> save stacked spec2d of all; %s' % coadd_name)

    return



# ---------------------------------------------------------------------------
# Image cutout helpers and 1-D spectral extraction
# ---------------------------------------------------------------------------

def _get_mosaic_cutout(
    coord: 'SkyCoord',
    band: str,
    size_arcsec: float,
    mosaic_dir: str,
    mosaic_filename_fmt: str,
    field: str,
) -> 'tuple[np.ndarray | None, WCS | None, float]':
    """
    Extract a 2D postage-stamp cutout from a large JWST mosaic FITS file on disk.

    Uses ``astropy.nddata.Cutout2D`` with memory-mapped I/O so only the
    required pixels are read from the (potentially very large) mosaic.

    Parameters
    ----------
    coord : SkyCoord
        Sky position of the target.
    band : str
        NIRCam filter name, e.g. ``'F444W'``.
    size_arcsec : float
        Side length of the square cutout in arcseconds.
    mosaic_dir : str
        Directory containing the mosaic FITS files.
    mosaic_filename_fmt : str
        ``str % (field, band.lower())`` format string for the mosaic filename.
    field : str
        Field identifier, e.g. ``'goods-s'`` or ``'goods-n'``.

    Returns
    -------
    data : 2D ndarray or None
        Cutout pixel values (NaN-padded at edges), or ``None`` on failure.
    cutout_wcs : WCS or None
        WCS of the returned cutout sub-image.
    pixscale : float
        Mosaic pixel scale in arcsec/pixel (``np.nan`` on failure).
    """
    import os as _os
    from astropy.wcs import WCS
    from astropy.nddata import Cutout2D

    fname = _os.path.join(mosaic_dir, mosaic_filename_fmt % (field, band.lower()))
    try:
        with fits.open(fname, memmap=True) as hdul:
            for hdu in hdul:
                if hdu.data is not None and hdu.data.ndim == 2:
                    wcs = WCS(hdu.header, naxis=2)
                    pixscale = float(np.abs(wcs.pixel_scale_matrix).max() * 3600)
                    size_pix = int(np.ceil(size_arcsec / pixscale))
                    cutout = Cutout2D(hdu.data, coord, size_pix, wcs=wcs,
                                     mode='partial', fill_value=np.nan)
                    return np.array(cutout.data, dtype=float), cutout.wcs, pixscale
    except Exception as exc:
        print('  [mosaic cutout] %s %s failed: %s' % (field, band, exc))
        return None, None, np.nan

    print('  [mosaic cutout] %s %s: no 2D data found' % (field, band))
    return None, None, np.nan


def _query_jades_cutout(
    coord: 'SkyCoord',
    band: str,
    size_arcsec: float,
    jades_field: str,
    jades_dr: str,
    version: str,
    target: str,
) -> 'tuple[np.ndarray | None, WCS | None, float]':
    """
    Fetch a 2D postage-stamp cutout from a JADES HLSP mosaic on MAST.

    Uses ``astrocut.FITSCutout`` with a MAST cloud URI to stream only the
    required pixels — no full-mosaic download needed.  This is an optional
    streaming alternative to :func:`_get_mosaic_cutout`; the 1D extraction
    worker currently uses :func:`_get_mosaic_cutout` instead.

    Parameters
    ----------
    coord : SkyCoord
        Sky position of the target.
    band : str
        NIRCam filter name, e.g. ``'F444W'``.
    size_arcsec : float
        Side length of the square cutout in arcseconds.
    jades_field : str
        JADES field name, e.g. ``'goods-s-deep'`` or ``'goods-n'``.
    jades_dr : str
        JADES data release tag, e.g. ``'dr1'``, ``'dr2'``, ``'dr3'``.
    version : str
        JADES mosaic version tag, e.g. ``'v1.0'``, ``'v2.0'``.
    target : str
        Target name of the field, e.g. ``'goods-s'`` or ``'goods-n'``.

    Returns
    -------
    data : 2D ndarray or None
        Cutout pixel values (NaN-padded at edges), or ``None`` on failure.
    cutout_wcs : WCS or None
        WCS of the returned cutout sub-image.
    pixscale : float
        Mosaic pixel scale in arcsec/pixel (``np.nan`` on failure).
    """
    from astrocut import FITSCutout
    from astropy.wcs import WCS
    import astropy.units as _u

    mast_uri = (
        'mast:HLSP/jades/{dr}/{target}/images/'
        'hlsp_jades_jwst_nircam_{field}_{band}_{version}_drz.fits'
    ).format(dr=jades_dr, target=target, field=jades_field, 
             band=band.lower(), version=version)

    try:
        cutout_obj = FITSCutout(
            input_files=[mast_uri],
            coordinates=coord,
            cutout_size=size_arcsec * _u.arcsec,
            extension=1,          # SCI extension only
            verbose=False,
        )
        hdul = cutout_obj.fits_cutouts[0]
    except Exception as exc:
        print('  [JADES cutout] %s %s failed: %s' % (jades_field, band, exc))
        return None, None, np.nan

    for hdu in hdul:
        if hdu.data is not None and hdu.data.ndim == 2:
            cutout_wcs = WCS(hdu.header, naxis=2)
            pixscale   = float(np.abs(cutout_wcs.pixel_scale_matrix).max() * 3600)
            return np.array(hdu.data, dtype=float), cutout_wcs, pixscale

    print('  [JADES cutout] %s %s: no 2D data found in result' % (jades_field, band))
    return None, None, np.nan


def extract_1d_spec_worker(
    spec2d_path: str,
    extraction_dir: str,
    do_boxcar: bool,
    grism_conf: "GrismConf",
    image_mosaic_dir: str = "",
    image_mosaic_filename_fmt: str = "hlsp_jades_jwst_nircam_%s_%s_v5.0_drz.fits",
    image_mosaic_rgb_bands: list = None,
    image_mosaic_field: str = 'goods-s',
    plot_dir: str = None,
) -> None:
    """
    Extract a 1D spectrum from a co-added 2D grism spectrum FITS file.

    Loads the 2D co-add, estimates an optimal spatial (cross-dispersion) profile
    from an image cutout read from a local mosaic FITS file on disk,
    performs optimal and/or boxcar extraction to produce a 1D spectrum, and
    saves a diagnostic plot (PDF) and the 1D spectrum table (FITS) alongside
    the input file.

    Parameters
    ----------
    spec2d_path : str
        Path to the co-added 2D spectrum FITS file, including all coadds.
    extraction_dir : str
        Directory to save the extracted 1D spectrum.
    do_boxcar : bool
        If True, skip optimal extraction and use boxcar only.
    grism_conf : GrismConf
        Grism configuration object containing dispersion and sensitivity
        calibration, used to convert pixel columns to wavelengths.
    image_mosaic_dir : str
        Directory containing the large JWST mosaic FITS files on disk.
    image_mosaic_filename_fmt : str
        ``str % (field, band.lower())`` format string for the mosaic filename.
    image_mosaic_rgb_bands : list of str
        Three-element list ``[blue_band, green_band, red_band]`` selecting the
        NIRCam filters for the B/G/R channels of the diagnostic thumbnail.
        Defaults to ``['F090W', 'F200W', 'F444W']``.
    image_mosaic_field : str
        Field identifier inserted into *image_mosaic_filename_fmt*
        (e.g. ``'goods-s'`` or ``'goods-n'``).
    plot_dir : str or None
        Directory for saving the diagnostic PDF plot.  If ``None``, the plot
        is saved alongside the input ``spec2d_path`` file.
    """
    if image_mosaic_rgb_bands is None:
        image_mosaic_rgb_bands = ['F090W', 'F200W', 'F444W']
    spec2d_fits = fits.open(spec2d_path)

    # --- Source metadata ---
    source_id = spec2d_fits[0].header['ID']
    filter_   = spec2d_fits[0].header['FILTER']
    obs_pa    = spec2d_fits[0].header['GS_V3_PA']  # V3 PA needed for morphology orientation

    if filter_ == 'F444W':
        wave_range = np.array([3.8, 5.1]);   mag_keyword = 'F444W_mag'
    elif filter_ == 'F322W2':
        wave_range = np.array([2.35, 4.1]);  mag_keyword = 'F356W_mag'
    elif filter_ == 'F356W':
        wave_range = np.array([3.05, 4.0]);  mag_keyword = 'F356W_mag'
    else:
        raise ValueError('filter %s not recognized' % filter_)

    source_mag = spec2d_fits[0].header[mag_keyword]
    try:
        source_mag = float(source_mag)
    except (TypeError, ValueError):
        source_mag = 99.
    if np.isnan(source_mag):
        source_mag = 99.
        spec2d_fits[0].header[mag_keyword] = 99.
    source_flux_mJy = 10 ** (-0.4 * source_mag) * 3631e3

    N_rows, N_cols = spec2d_fits[0].header['N_R'], spec2d_fits[0].header['N_C']
    if N_cols > N_rows:
        obs_pa -= 90

    tb_stats = Table(spec2d_fits['STATS'].data)
    N_A = np.sum(tb_stats['module'] == 'A')
    N_B = np.sum(tb_stats['module'] == 'B')

    RA, DEC = spec2d_fits[0].header['RA0'], spec2d_fits[0].header['DEC0']
    coord   = SkyCoord(RA, DEC, unit=(u.deg, u.deg))

    # --- Spatial profile estimation from image cutout ---
    # NOTE: A and B are in arcsec (JADES DR5), PA in degrees (N→E); we rotate PA
    # to be relative to the dispersion direction using obs_pa (V3 PA).
    profile_name    = 'none'
    profile_spline  = None
    aper_corr       = 1.0
    chisq_full      = np.nan
    pixscale_profile = np.nan
    try:
        A, B, PA = spec2d_fits[0].header['A'], spec2d_fits[0].header['B'], spec2d_fits[0].header['PA']
        A, B = np.clip(A, 0.06, 10), np.clip(B, 0.06, 10)
        PA   = (PA + obs_pa) % 180

        # Read the grism-filter mosaic from disk for the spatial profile
        size_arcsec_profile = NIRCAM_LW_PIXSCALE * spec2d_fits['spec2d'].data.shape[0] * 1.5
        DI_cutout, _, pixscale_profile = _get_mosaic_cutout(
            coord, filter_, size_arcsec_profile,
            image_mosaic_dir, image_mosaic_filename_fmt, image_mosaic_field,
        )
        if DI_cutout is None:
            raise KeyError('JADES cutout unavailable')

        DI_cutout_rot = ndimage.rotate(np.nan_to_num(DI_cutout), obs_pa, reshape=False)
        DI_cutout_rot[DI_cutout_rot == 0] = np.nan
        profile_data  = np.nansum(DI_cutout_rot, axis=1)
        profile_name  = 'image_collapse'

        # Fit Gaussian to the collapsed spatial profile
        try:
            popt, _ = optimize.curve_fit(
                f=gauss,
                xdata=np.arange(len(profile_data)),
                ydata=profile_data / np.max(profile_data),
                p0=[len(profile_data) // 2, 1, 1],
            )
            gauss_fit    = gauss(np.arange(len(profile_data)), *popt)
            profile_norm = profile_data / np.max(profile_data)
            chisq_full   = np.sum((gauss_fit - profile_norm) ** 2)
            idx_center   = np.arange(int(len(profile_data) * 0.33), int(len(profile_data) * 0.66))
            chisq_center = np.sum(
                (gauss_fit[idx_center] - profile_norm[idx_center] / np.max(profile_norm[idx_center])) ** 2
            )
        except (ValueError, RuntimeError):
            chisq_full, chisq_center = np.nan, np.nan
            popt = np.array([np.nan, np.nan, np.nan])

        center_ok = np.abs(popt[0] - len(profile_data) / 2) <= len(profile_data) / 10

        if (chisq_full ** 0.5 > 0.5) and (chisq_center ** 0.5 < 0.5) and center_ok:
            # Full fit poor but center well-described: use best-fit Gaussian
            profile_data = gauss(np.arange(len(profile_data)), *popt)
            profile_name = 'image_collapse_gaussian'
        elif (np.sum(np.isfinite(profile_data) & (profile_data != 0)) == 0) or \
                np.isnan(chisq_center) or (chisq_full ** 0.5 > 0.5) or not center_ok:
            # Image model unavailable or unreliable: fall back to parametric Sersic
            ny, nx = DI_cutout.shape
            xx, yy = np.meshgrid(np.arange(nx), np.arange(ny))
            sersic = models.Sersic2D(
                x_0=nx / 2, y_0=ny / 2,
                amplitude=1., r_eff=A / pixscale_profile, n=1,
                ellip=1 - B / A, theta=np.deg2rad(PA),
            )
            profile_data = np.sum(sersic(xx, yy), axis=1)
            profile_name = 'image_SExtractor_model'

        # Build interpolated spline mapping profile pixel coords onto spectrum pixel coords
        n_spec_rows = spec2d_fits['spec2d'].data.shape[0]
        profile_x   = (np.arange(len(profile_data)) - len(profile_data) // 2) * pixscale_profile / NIRCAM_LW_PIXSCALE + n_spec_rows // 2
        try:
            profile_spline = interpolate.UnivariateSpline(
                profile_x, profile_data / np.max(profile_data), s=0, k=1, ext='zeros'
            )
        except Exception:
            profile_spline = interpolate.UnivariateSpline(
                np.arange(31), np.concatenate((np.zeros(13), np.ones(5), np.zeros(13))),
                s=0, k=1, ext='zeros',
            )
            profile_name = 'boxcar_5pix'

    except KeyError:
        profile_name = 'none'
        do_boxcar    = True

    # --- Continuum contamination flag ---
    try:
        cont_mag_limit = np.log10(9 / 10 * 2 * (spec2d_fits[0].header['EFFEXPTM'] / 1e4) ** -0.5 / 3631e6) * -2.5
    except KeyError:
        cont_mag_limit = 23.0
    is_cont = source_mag < cont_mag_limit

    # --- Load 2D spectral arrays ---
    bunit = spec2d_fits['spec2d'].header.get('BUNIT')
    spec2d = spec2d_fits['spec2d'].data
    line2d = spec2d_fits['line2d'].data
    wht2d  = spec2d_fits['wht2d'].data
    wave   = spec2d_fits[1].header['WAVE_1'] + \
             np.arange(spec2d_fits[1].header['NAXIS1']) * spec2d_fits[1].header['D_WAVE']

    # --- Re-estimate continuum via median filtering; keep if it reduces noise ---
    line2d_orig   = line2d.copy()
    highSN_mask   = line2d_orig / wht2d ** -0.5 > 2.0
    if bunit == "mJy":
        spec2d_counts = spec2d * (grism_conf.get_sensitivity("A", "R")(wave) * 1e-3)
    elif bunit == "DN/s":
        spec2d_counts = spec2d.copy()
    else:
        raise ValueError('Unrecognized BUNIT %s in spec2d header' % bunit)
    valid_cols    = np.where(np.sum(np.isnan(spec2d_counts), axis=0) != len(spec2d_counts))[0]

    spec2d_counts_medflt = spec2d_counts.copy()
    spec2d_counts_medflt[:, valid_cols] = ndimage.median_filter(
        spec2d_counts[:, valid_cols], footprint=np.ones((1, 150)), mode='reflect'
    )
    spec2d_counts[highSN_mask] = spec2d_counts_medflt[highSN_mask]
    spec2d_counts_medflt_narrow = spec2d_counts_medflt.copy()
    spec2d_counts_medflt_narrow[:, valid_cols] = ndimage.median_filter(
        spec2d_counts[:, valid_cols], footprint=np.ones((1, 50)), mode='reflect'
    )
    if bunit == "mJy":
        line2d_new = spec2d - np.nan_to_num(spec2d_counts_medflt_narrow / (grism_conf.get_sensitivity("A", "R")(wave) * 1e-3))
    else:
        line2d_new = spec2d - np.nan_to_num(spec2d_counts_medflt_narrow)
    if sigma_clipped_stats(line2d_new, sigma=2)[2] <= sigma_clipped_stats(line2d_orig, sigma=2)[2]:
        line2d = line2d_new

    if np.sum(~np.isnan(np.nansum(spec2d, axis=0))) < 200:
        return

    # --- Extraction aperture (boxcar half-width = 2 pixels) ---
    yc   = spec2d.shape[0] // 2
    aper = 2
    sl   = slice(yc - aper, yc + aper + 1)

    # --- Optimal extraction ---
    wave_opt = spec1d_opt = spec1d_cont_opt = unc1d_opt = None
    if not do_boxcar:
        profile_1d  = profile_spline(np.arange(line2d.shape[0], dtype=float))
        profile_1d /= np.nansum(profile_1d)
        denom           = np.nansum(wht2d.T * profile_1d ** 2, axis=1)
        spec1d_cont_opt = np.nansum((spec2d * wht2d).T * profile_1d, axis=1) / denom
        spec1d_opt      = np.nansum((line2d * wht2d).T * profile_1d, axis=1) / denom
        unc1d_opt       = denom ** -0.5
        good_cols_opt   = np.where(np.sum(wht2d == 0, axis=0) == 0)[0]
        wave_opt        = wave[good_cols_opt]
        spec1d_cont_opt = spec1d_cont_opt[good_cols_opt]
        spec1d_opt      = spec1d_opt[good_cols_opt]
        unc1d_opt       = unc1d_opt[good_cols_opt]

    # --- Boxcar extraction with aperture correction ---
    spec1d_cont_box = np.nansum(spec2d[sl], axis=0)
    spec1d_box      = np.nansum(line2d[sl], axis=0)
    unc1d_box       = np.nansum(wht2d[sl] ** -1, axis=0) ** 0.5

    if profile_spline is not None:
        profile_full = profile_spline(np.arange(line2d.shape[0], dtype=float))
        profile_full = profile_full / np.nansum(profile_full)
        aper_corr    = np.sum(profile_full) / np.sum(profile_full[sl])
    spec1d_cont_box *= aper_corr
    spec1d_box      *= aper_corr
    unc1d_box       *= aper_corr

    good_cols_box   = np.where(np.sum(wht2d[sl] == 0, axis=0) == 0)[0]
    wave_box        = wave[good_cols_box]
    spec1d_cont_box = spec1d_cont_box[good_cols_box]
    spec1d_box      = spec1d_box[good_cols_box]
    unc1d_box       = unc1d_box[good_cols_box]

    # --- Select extraction method ---
    if do_boxcar:
        wave, spec1d_cont, spec1d, unc1d = wave_box, spec1d_cont_box, spec1d_box, unc1d_box
    else:
        snr_opt = np.nanmedian((spec1d_cont_opt / unc1d_opt)[unc1d_opt < np.nanmedian(unc1d_opt) * 2])
        snr_box = np.nanmedian((spec1d_cont_box / unc1d_box)[unc1d_box < np.nanmedian(unc1d_box) * 2])
        if (snr_opt > snr_box) and is_cont:
            wave, spec1d_cont, spec1d, unc1d = wave_opt, spec1d_cont_opt, spec1d_opt, unc1d_opt
        else:
            do_boxcar = True
            wave, spec1d_cont, spec1d, unc1d = wave_box, spec1d_cont_box, spec1d_box, unc1d_box

    if np.nanmedian(spec1d[spec1d != 0]) < 0:
        spec1d -= np.nanmedian(spec1d[spec1d != 0])

    # =========================================================================
    # Diagnostic plot
    # =========================================================================
    plt.close()
    e, b  = 0.1, 0.9   # edge and bottom margins [inches]
    W, H  = 16 + e * 3, 7 + e * 3 + b
    fig   = plt.figure(figsize=(W / H * 6, 6))
    ax_im = fig.add_axes([e/W,         (5 + b + 2*e)/H, 2/W,  2/H])  # direct image
    ax_2d = fig.add_axes([(2 + e*2)/W, (5 + b + 2*e)/H, 14/W, 2/H])  # 2D spectrum (full)
    ax_li = fig.add_axes([(2 + e*2)/W, (3 + b +   e)/H, 14/W, 2/H])  # 2D spectrum (line)
    ax_1d = fig.add_axes([(2 + e*2)/W, b/H,             14/W, 3/H])  # 1D spectrum
    ax    = [ax_2d, ax_li, ax_1d]

    D_WAVE   = spec2d_fits[1].header['D_WAVE']
    WAVE_1   = spec2d_fits[1].header['WAVE_1']
    aspect   = (np.diff(wave_range)[0] - 0.1) / D_WAVE / (spec2d.shape[0] - 1) / 7
    xticks   = np.arange(np.ceil((wave_range[0] + 0.05) * 10) / 10, wave_range[1] - 0.04, 0.1)
    xtick_px = (xticks - WAVE_1) / D_WAVE
    xlim_px  = ((wave_range[0] + 0.05 - WAVE_1) / D_WAVE, (wave_range[1] - WAVE_1 - 0.05) / D_WAVE)

    # 2D full spectrum
    vmin_2d, vmax_2d = ZScaleInterval().get_limits(spec2d[:, 100:-100])
    ax[0].imshow(spec2d, aspect=aspect, vmin=vmin_2d, vmax=vmax_2d,
                 cmap=plt.cm.gist_gray_r, origin='lower')
    ax[0].set(ylim=(0.5, spec2d.shape[0] - 0.5), aspect=aspect, xticks=[], xlim=xlim_px)
    ax[0].set_yticks([spec2d.shape[0] / 2.]); ax[0].set_yticklabels([''])
    ax[0].set_xticks(xtick_px); ax[0].set_xticklabels([])

    # 2D continuum-subtracted spectrum
    vmin_li, vmax_li = ZScaleInterval().get_limits(line2d[:, 100:-100])
    ax[1].imshow(line2d, vmin=vmin_li / 2., vmax=vmax_li,
                 cmap=plt.cm.gist_gray_r, origin='lower')
    ax[1].set(ylim=(0.5, spec2d.shape[0] - 0.5), aspect=aspect, xticks=[], xlim=xlim_px)
    ax[1].set_yticks([spec2d.shape[0] / 2.]); ax[1].set_yticklabels([''])
    ax[1].set_xticks(xtick_px); ax[1].set_xticklabels([])

    for ax_ in ax[:2]:
        ax_.yaxis.set_tick_params(width=1.5, size=4, right=True)
    if is_cont:
        ax[0].axhline(yc + aper + 1.5, color='w', ls='--', dashes=(4, 4))
        ax[0].axhline(yc - aper - 0.5, color='w', ls='--', dashes=(4, 4))
    ax[1].axhline(yc + aper + 1.5, color='w', ls='--', dashes=(4, 4))
    ax[1].axhline(yc - aper - 0.5, color='w', ls='--', dashes=(4, 4))

    # 1D spectrum
    kw1d = dict(lw=1.5, drawstyle='steps-mid')
    ax[2].plot(wave, ndimage.gaussian_filter1d(spec1d, 0.6), color='k', zorder=100, **kw1d)
    if is_cont:
        ax[2].plot(wave, ndimage.gaussian_filter1d(spec1d_cont, 0.6), color='dimgrey', **kw1d)
        ymax_1d = np.nanpercentile(spec1d_cont[np.isfinite(spec1d_cont) & (spec1d_cont != 0)], 95) * 1.25
    else:
        ymax_1d = np.nanpercentile(spec1d[np.isfinite(spec1d) & (spec1d != 0)], 95) * 1.5
    ax[2].axhline(0, color='grey', ls='--')
    ax[2].set(xlim=(wave_range[0] + 0.05, wave_range[1] - 0.05), xticks=xticks,
              xlabel='Observed Wavelength (µm)', ylabel='Flux Density [{}]'.format(bunit))
    ax[2].set_ylim(np.clip(vmin_li * 2.0, -0.035, 0), np.clip(ymax_1d, 0.015, 1e8))

    # Annotations
    corner_text(ax[0], loc=2, s='ID%s' % source_id, weight='semibold', color='r', fontsize=20, edge=5e-3)
    corner_text(ax[0], loc=1, s='%s=%.2f' % (mag_keyword, source_mag), color='r', fontsize=15, edge=5e-3)
    if N_A == 0:
        corner_text(ax[0], s='modB', loc=3, color='r', fontsize=14, edge=5e-3,
                    path_effects=[pe.withStroke(linewidth=2.5, foreground='w')])
    elif N_B == 0:
        corner_text(ax[0], s='modA', loc=3, color='r', fontsize=14, edge=5e-3,
                    path_effects=[pe.withStroke(linewidth=2.5, foreground='w')])
    else:
        corner_text(ax[0], s='modA:%d / modB:%d' % (N_A, N_B), loc=3, color='r', fontsize=14, edge=5e-3,
                    path_effects=[pe.withStroke(linewidth=2.5, foreground='w')])
    corner_text(ax[1], loc=3, s='Continuum Subtracted', color='r', fontsize=15, edge=5e-3,
                path_effects=[pe.withStroke(linewidth=2.5, foreground='w')])
    corner_text(ax[1], loc=4, s='(%.5f, %.5f)' % (RA, DEC), color='r', fontsize=15, edge=5e-3)
    corner_text(ax[2], loc=1, s='$f_\\mathrm{%s}$=%.3f mJy' % (mag_keyword.split('_')[0], source_flux_mJy),
                fontsize=15, color='r', edge=5e-3, path_effects=[pe.withStroke(linewidth=2.5, foreground='w')])

    # Redshift and emission-line markers
    if 'z_grism' in spec2d_fits[0].header:
        z = spec2d_fits[0].header['z_grism']
        corner_text(ax[2], loc=4, s='z=%.3f' % z, color='r', fontsize=15, edge=5e-3, zorder=999,
                    path_effects=[pe.withStroke(linewidth=2.5, foreground='w')])
        kw_vline = dict(ymin=0., ymax=1., zorder=-5, lw=5, alpha=0.5, color='skyblue')
        line_names = np.array([
            r'[O$\,$II]',    r'H$\rm\beta$',   r'[O$\,$III]',   r'[O$\,$III]',  r'H$\rm\alpha$',
            r'[N$\,$II]',    r'[S$\,$II]',      r'[S$\,$III]',   r'[S$\,$III]',
            r'Pa$\rm\delta$', r'He$\,$I',       r'Pa$\rm\gamma$', r'[Fe$\,$II]', r'Pa$\rm\beta$',
            r'[Fe$\,$II]',   r'Pa$\rm\alpha$',  r'He$\,$I',      'H$_2$',        r'Br$\rm\gamma$',
            'H$_2$',         'H$_2$',           r'Br$\rm\beta$', 'H$_2$',        'PAH',
            r'Pf$\,$8',      r'Br$\rm\alpha$',
        ])
        line_waves = np.array([
             3728.5,  4862.67,  4960.295,  5008.24,  6564.61,  6585.27,  6725.48,
             9071.1,  9533.21, 10052.1,   10833.3,  10941.0,  12570.2,  12821.5,
            16440.5, 18756.0,  20592.5,   21223.8,  21661.0,  24072.6,  24243.6,
            26258.4, 28032.6,    32900,   37405.2,  40522.3,
        ])
        w0, w1 = wave_range[0] + 0.05, wave_range[1] - 0.05
        for name, lam in zip(line_names[::-1], line_waves[::-1]):
            w_obs = (1 + z) * lam / 1e4
            if not (w0 < w_obs < w1):
                continue
            ax[2].axvline(w_obs, **kw_vline)
            ax[2].text(w_obs, ax[2].get_ylim()[0] * 0.25 + ax[2].get_ylim()[1] * 0.75,
                       s=name, zorder=-2, color='b', ha='center', va='center',
                       rotation=90, fontsize=14,
                       path_effects=[pe.withStroke(linewidth=5, foreground='w')])
    elif 'z_a' in spec2d_fits[0].header:
        corner_text(ax[2], loc=4,
                    s='z=%.2f (%.2f\u2013%.2f)' % (spec2d_fits[0].header['z_a'],
                                                    spec2d_fits[0].header.get('z_16', float('nan')),
                                                    spec2d_fits[0].header.get('z_84', float('nan'))),
                    color='r', fontsize=15, edge=5e-3, zorder=999,
                    path_effects=[pe.withStroke(linewidth=2.5, foreground='w')])

    # Direct image cutout panel — RGB composite from mosaic files on disk
    # image_mosaic_rgb_bands = [blue_band, green_band, red/LW_band]
    band_blue, band_green, band_lw = image_mosaic_rgb_bands
    size_arcsec_di = NIRCAM_LW_PIXSCALE * spec2d.shape[0] * 1.5
    DI_LW,    _, pixscale_lw = _get_mosaic_cutout(coord, band_lw,    size_arcsec_di, image_mosaic_dir, image_mosaic_filename_fmt, image_mosaic_field)
    DI_green, _, _           = _get_mosaic_cutout(coord, band_green, size_arcsec_di, image_mosaic_dir, image_mosaic_filename_fmt, image_mosaic_field)
    DI_blue,  _, _           = _get_mosaic_cutout(coord, band_blue,  size_arcsec_di, image_mosaic_dir, image_mosaic_filename_fmt, image_mosaic_field)

    if DI_LW is None:
        DI_LW = np.zeros((10, 10)); pixscale_lw = NIRCAM_LW_PIXSCALE
    hf_box  = int(NIRCAM_LW_PIXSCALE * (spec2d.shape[0] - 1) / 2. / pixscale_lw)
    di_size = DI_LW.shape[0]

    DI_LW_rot = ndimage.rotate(np.nan_to_num(DI_LW), obs_pa, reshape=False)
    DI_LW_rot[DI_LW_rot == 0] = np.nan

    if DI_green is not None and DI_blue is not None:
        DI_green_rot = ndimage.rotate(np.nan_to_num(DI_green), obs_pa, reshape=False)
        DI_blue_rot  = ndimage.rotate(np.nan_to_num(DI_blue),  obs_pa, reshape=False)
        DI_green_rot[DI_green_rot == 0] = np.nan
        DI_blue_rot[DI_blue_rot == 0]   = np.nan
        DI_rgb = np.dstack((DI_LW_rot, DI_green_rot, DI_blue_rot))
    else:
        DI_rgb = np.dstack((DI_LW_rot, np.zeros_like(DI_LW_rot), np.zeros_like(DI_LW_rot)))
    DI_rgb = np.clip((np.log10(DI_rgb / np.clip(np.nanpercentile(DI_rgb, 99.5), 1, 10) + 1e-2) + 2.0) / 2.0, 0, 1)

    if (N_A == 0) and (N_cols == 0):
        DI_LW_rot = DI_LW_rot[:, ::-1]; DI_rgb = DI_rgb[:, ::-1, :]
    if N_rows == 0:
        DI_LW_rot = DI_LW_rot[:, ::-1]; DI_rgb = DI_rgb[:, ::-1, :]

    try:
        di_vmin, di_vmax = ZScaleInterval().get_limits(DI_LW_rot[np.isfinite(DI_LW_rot)])
        di_vmax *= 2
    except IndexError:
        di_vmin, di_vmax = 0.0, 0.1

    ax_im.imshow(np.nan_to_num(DI_rgb), cmap=plt.cm.gist_heat, origin='lower')
    ax_im.set(xlim=(di_size // 2 - hf_box - 1, di_size // 2 + hf_box - 1),
              ylim=(di_size // 2 - hf_box - 1, di_size // 2 + hf_box - 1))
    ax_im.set_xticks([]); ax_im.set_yticks([])
    band_label = '-'.join(b[1:4] for b in [band_blue, band_green, band_lw])
    corner_text(ax_im, loc=2, s=band_label, color='w', fontsize=13, weight='semibold')

    print(source_id, '%s=%.2f' % (mag_keyword, source_mag),
          'chi_img_model=%.3f' % chisq_full ** 0.5, profile_name, 'do_boxcar=', do_boxcar)

    # =========================================================================
    # Save outputs
    # =========================================================================
    fig_output_path = os.path.join(plot_dir, 
                 os.path.basename(spec2d_path).replace('.fits', '.pdf'))
    plt.savefig(fig_output_path, dpi=150)

    tb_box = Table(
        names=['wavelength_um', 'box_spec1d_mJy', 'box_line1d_mJy', 'box_fluxerr_mJy'],
        data=[wave_box, spec1d_cont_box, spec1d_box, unc1d_box],
    )
    if profile_name != 'none' and wave_opt is not None:
        tb_opt = Table(
            names=['wavelength_um', 'opt_spec1d_mJy', 'opt_line1d_mJy', 'opt_fluxerr_mJy'],
            data=[wave_opt, spec1d_cont_opt, spec1d_opt, unc1d_opt],
        )
        tb_1d = join(tb_opt, tb_box, keys='wavelength_um', join_type='outer')
    else:
        tb_1d = tb_box

    for col in tb_1d.colnames:
        tb_1d[col].fill_value = np.nan
        tb_1d[col].info.format = '.5f' if 'mJy' in col else '.4f'
    tb_1d = tb_1d.filled()

    tb_1d.meta['comments'] = [
        '-' * 70,
        '1D spectrum extracted at y_c=%.1f with aperture height = %.1f pix' % (yc, aper * 2 + 1),
        ('I only subtracted common grism sky background. Contaminants are not subtracted.'
         if is_cont else
         'Extracted from 2D grism images that have been continuum/background-subtracted.'),
        'Be careful about potential contaminant & aperture loss.',
        'Produced by F. Sun (CfA | Harvard & Smithsonian, %s)' % time.strftime('%Y/%m/%d', time.localtime()),
        '-' * 70,
    ]
    out_path = os.path.join(extraction_dir,
                os.path.basename(spec2d_path).replace('spec_2d_', 'spec_1d_'))
    fits_1d  = fits.BinTableHDU(tb_1d)
    for card in spec2d_fits[0].header.cards[4:]:
        fits_1d.header['HIERARCH ' + card[0]] = (card[1], card[2])
    fits_1d.header['N_A']      = (N_A,          'number of mod A exposures')
    fits_1d.header['N_B']      = (N_B,          'number of mod B exposures')
    fits_1d.header['boxcar']   = (do_boxcar,     'extracted using boxcar method')
    fits_1d.header['y_c']      = (yc,            'y_center [pix] of aperture in 2D spectrum')
    fits_1d.header['aper']     = (aper * 2 + 1,  'total aperture height [pix]')
    fits_1d.header['profile']  = (profile_name,  'image profile used for aperture correction')
    fits_1d.header['apercorr'] = (aper_corr,     'aperture correction factor')
    fits_1d.writeto(out_path, overwrite=True)
# ---------------------------------------------------------------------------
# 2-D emission-line cutout extraction
# ---------------------------------------------------------------------------

def store_all_2d_emline(
    cutout_list: list,
    xs_list: list,
    ys_list: list,
    theta_list: list,
    modules: list,
    pupils: list,
    paths: list,
    effexptm_list: list,
    gs_v3pa_list: list,
    output: str,
    coord: "SkyCoord",
    source_item: "Table.Row",
    grism_filter: str,
    wave_line: float,
    name_line: str,
    overwrite: bool = True,
) -> fits.HDUList:
    """
    Combine per-frame 2-D emission-line cutouts into a single FITS file.

    Each entry in ``cutout_list`` is a tuple ``(sci, line, wht, dq)`` of 2-D
    arrays (NxN, in DN/s).  One set of extensions is written per frame,
    following a layout similar to :func:`store_all_2d_spec`.

    Parameters
    ----------
    cutout_list:
        List of ``(sci, line, wht, dq)`` tuples for each frame.
    xs_list, ys_list:
        Emission-line pixel position in the native grism frame.
    theta_list:
        Local dispersion angle (rad) in the sky-aligned output frame.
    modules, pupils, paths:
        Detector module, pupil, and file path for each frame.
    effexptm_list, gs_v3pa_list:
        Exposure time and V3 position angle for each frame.
    output:
        Output file path.
    coord:
        Sky coordinate of the source.
    source_item:
        Single-row ``Table.Row`` whose columns are written to the primary header.
    grism_filter, wave_line, name_line:
        Filter, observed wavelength (µm), and name of the targeted emission line.
    overwrite:
        Overwrite the output file if it exists.

    Returns
    -------
    hdul : fits.HDUList
        The assembled ``HDUList`` (also written to disk when ``overwrite=True``).
    """
    hdul = fits.HDUList([fits.PrimaryHDU()])
    hdul[0].header["RA0"]      = (float(coord.ra.deg),   "Source RA (deg)")
    hdul[0].header["DEC0"]     = (float(coord.dec.deg),  "Source DEC (deg)")
    hdul[0].header["FILTER"]   = (grism_filter,          "Filter name")
    hdul[0].header["LINENAME"] = (name_line,             "Target emission line")
    hdul[0].header["LINEWAVE"] = (wave_line,             "Observed wavelength of line (um)")
    hdul[0].header["N_COADD"] = (len(cutout_list),       "Total number of frames")
    hdul[0].header["AUTHOR"]   = ("Jiachuan Xu",         "Author")
    hdul[0].header["TIME"]     = (
        time.strftime("%Y/%m/%d %H:%M:%S", time.localtime()), "Creation time"
    )
    if source_item is not None:
        hdul[0].header["COMMENTS"] = "Source catalog information below:"
        for col in source_item.colnames:
            val = source_item[col]
            if isinstance(val, np.ma.core.MaskedConstant):
                continue
            try:
                if isinstance(val, (str, np.str_)):
                    hdul[0].header["HIERARCH " + col] = str(val)
                else:
                    fval = float(val)
                    if np.isnan(fval):
                        hdul[0].header["HIERARCH " + col] = "nan"
                    elif np.isinf(fval):
                        hdul[0].header["HIERARCH " + col] = "inf"
                    else:
                        hdul[0].header["HIERARCH " + col] = val
            except Exception:
                pass

    for i, (sci, line, wht, dq) in enumerate(cutout_list):
        h_sci = fits.ImageHDU(sci.astype(np.float32), name="SCI-%d" % i)
        h_sci.header["XS"]       = (float(xs_list[i]),    "X position in grism frame (pix)")
        h_sci.header["YS"]       = (float(ys_list[i]),    "Y position in grism frame (pix)")
        h_sci.header["DISPANG"]    = (float(np.rad2deg(theta_list[i])), 
                                      "Dispersion angle w.r.t the cutout X axis (deg, CCW from +X)")
        h_sci.header["MODULE"]   = modules[i]
        h_sci.header["PUPIL"]    = pupils[i]
        h_sci.header["FRAME_PA"] = (float(0.0), 
                                    "Position angle of the cutout X axis (northward from west)")
        h_sci.header["DATAPATH"] = os.path.basename(paths[i])
        h_sci.header["EFFEXPTM"] = (float(effexptm_list[i]), "Effective exposure time (s)")
        h_sci.header["GS_V3_PA"] = (float(gs_v3pa_list[i]),  "V3 position angle (deg)")
        h_sci.header["BUNIT"]    = ("DN/s", "Brightness unit")
        hdul.append(h_sci)
        hdul.append(fits.ImageHDU(wht.astype(np.float32),    name="WHT-%d" % i))
        hdul.append(fits.ImageHDU(dq.astype(np.int16),       name="DQ-%d"  % i))
        hdul.append(fits.ImageHDU(line.astype(np.float32),   name="LINE-%d" % i))

    stats = Table(
        names=["index", "xs", "ys", "theta_deg", "module", "pupil",
               "EFFEXPTM", "GS_V3_PA", "datapath"],
        data=[
            list(range(len(cutout_list))),
            [float(v) for v in xs_list],
            [float(v) for v in ys_list],
            [float(np.rad2deg(t)) for t in theta_list],
            modules, pupils,
            [float(v) for v in effexptm_list],
            [float(v) for v in gs_v3pa_list],
            [os.path.basename(p) for p in paths],
        ],
    )
    hdul.append(fits.BinTableHDU(stats, name="STATS"))

    if overwrite:
        hdul.writeto(output, overwrite=True)
    return hdul


def _compute_grism_psf_frame(
    primary_hd: fits.Header,
    xs: float,
    ys: float,
    filter_name: str,
    oversample: int = 4,
    fov_pixels: int = 51,
) -> tuple[np.ndarray, int, str]:
    """
    Compute a super-sampled NIRCam PSF for one grism exposure frame.

    The PSF is treated as achromatic within the narrow spectral window of the
    emission-line cutout, so we compute an imaging PSF (CLEAR pupil) at the
    same filter and detector position as the grism exposure.  This is an
    excellent approximation because the grism element does not materially
    change the spatial PSF shape.

    In-flight wavefront sensing (WSS) OPD data are loaded from the MAST
    archive via ``stpsf.load_wss_opd_by_date`` when the observation date
    keyword ``DATE-BEG`` is present in *primary_hd*.  If that call fails
    (no internet, OPD not available, stpsf not installed), the function
    falls back gracefully to the default stpsf OPD or to a simple Gaussian
    approximation.

    Parameters
    ----------
    primary_hd:
        Primary FITS header of the grism level-1.5 file.  Must contain at
        least ``MODULE`` (``'A'`` or ``'B'``).  ``DETECTOR``, ``DATE-BEG``,
        and ``GS_V3_PA`` are used when available.
    xs, ys:
        Predicted pixel position of the emission line in the grism detector
        frame (0-indexed).
    filter_name:
        NIRCam LW filter, e.g. ``'F444W'``.
    oversample:
        Super-sampling factor relative to the native pixel scale.  The
        returned array has ``fov_pixels * oversample`` pixels on each side.
    fov_pixels:
        Field of view in *native* detector pixels.  Should equal the
        ``cutout_size`` used for the science cutout.

    Returns
    -------
    psf_arr : ndarray, shape (fov_pixels*oversample, fov_pixels*oversample)
        Normalised PSF (sums to 1.0) at the oversampled pixel scale.
    oversample_used : int
        Actual oversample factor used (equals *oversample* on success,
        1 in the Gaussian fallback).
    stpsf_version : str
        ``stpsf.__version__`` string, or ``'gaussian_fallback'``.
    """
    module   = primary_hd.get("MODULE", "A")
    det_name = primary_hd.get("DET_NAME", "NRC%s5" % module)
    date_obs = primary_hd.get("DATE-BEG", primary_hd.get("DATE-OBS", None))
    assert (xs >= VALID_X_LEFT) and (xs <= VALID_X_RIGHT) and (ys >= VALID_Y_BOTTOM) and (ys <= VALID_Y_TOP), \
         "Source position (%.1f, %.1f) is too close to the detector edge for reliable PSF computation" % (xs, ys)
    try:
        import stpsf
        stpsf_ver = stpsf.__version__
    except ImportError:
        stpsf_ver = None

    if stpsf_ver is not None:
        try:
            nc = stpsf.NIRCam()
            nc.filter          = filter_name
            nc.detector        = det_name
            nc.detector_position = (float(xs), float(ys))
            # Use imaging (CLEAR) pupil — grism element does not change spatial PSF
            nc.pupil_mask      = None

            # In-flight OPD: try to load WSS measurement for this date
            if date_obs is not None:
                try:
                    nc.load_wss_opd_by_date(date_obs, plot=False, verbose=False)
                except Exception as _e_opd:
                    pass  # fall back to default OPD already set on nc

            psf_hdul  = nc.calc_psf(
                #monochromatic=?,
                oversample=oversample,
                fov_pixels=fov_pixels,
                normalize="last",
            )
            psf_arr   = psf_hdul["OVERSAMP"].data.copy()
            psf_arr  /= psf_arr.sum()
            return psf_arr, oversample, stpsf_ver

        except Exception as _e_stpsf:
            pass  # fall through to Gaussian fallback

    # --- Gaussian fallback -------------------------------------------------------
    # Approximate NIRCam LW PSF FWHM: λ/D ≈ 1.22 * (wave_um / 6.5 m) in arcsec,
    # converted to native pixels.  Use the filter central wavelength as a proxy.
    _FILTER_WAVE = {
        "F277W": 2.77, "F356W": 3.56, "F410M": 4.10,
        "F444W": 4.44, "F460M": 4.60, "F470N": 4.71,
    }
    wave_um = _FILTER_WAVE.get(filter_name, 4.0)
    # Diffraction-limited FWHM in arcsec, converted to native pix
    fwhm_arcsec = 1.22 * wave_um * 1e-6 / 6.5 * (180 / np.pi) * 3600
    fwhm_pix    = fwhm_arcsec / (NIRCAM_LW_PIXSCALE * oversample) # supersampled pixels
    sigma_over  = fwhm_pix / (2.0 * np.sqrt(2.0 * np.log(2.0)))   # supersampled σ

    N_over = fov_pixels * oversample  # Gaussian fallback is at super-sampled scale
    y_g, x_g = np.mgrid[0:N_over, 0:N_over]
    c = (N_over - 1) / 2.0
    psf_arr = np.exp(-((x_g - c) ** 2 + (y_g - c) ** 2) / (2.0 * sigma_over ** 2))
    psf_arr /= psf_arr.sum()
    return psf_arr, oversample, "gaussian_fallback"


# ---------------------------------------------------------------------------
# Helpers for extract_2d_emline_worker_drizzle
# ---------------------------------------------------------------------------

def _make_sky_wcs(ra, dec, pixscale_arcsec, n_pix):
    """Sky-aligned TAN WCS centred at (ra, dec); +X = west, +Y = north."""
    wcs_header = fits.Header(
        {"NAXIS": 2, 
        "NAXIS1": n_pix, 
        "NAXIS2": n_pix,
        "WCSAXES": 2,
        "CRPIX1": n_pix / 2.0 + 1.0,
        "CRPIX2": n_pix / 2.0 + 1.0,
        "CRVAL1": ra,
        "CRVAL2": dec,
        "CTYPE1": 'RA---TAN-SIP',
        "CTYPE2": 'DEC--TAN-SIP',
        "CUNIT1": 'deg',
        "CUNIT2": 'deg',
        "CD1_1": -pixscale_arcsec / 3600.0,
        "CD1_2": 0.0,
        "CD2_1": 0.0,
        "CD2_2": pixscale_arcsec / 3600.0}
    )
    wcs = WCS(wcs_header)
    return wcs, wcs.to_header()


def _build_cutout_and_pixmap(image, emline, weight, data_quality,
                              x0, y0, xs, ys, cutout_size, wcs_di, out_wcs):
    """
    Take a padded cutout centred at (xs, ys) and build the drizzle pixel map.

    The cutout is padded to ``ceil(cutout_size * sqrt(2) / 2) + 1`` half-width
    so that the sky-aligned output square is always fully covered at any roll angle.

    Returns a dict with keys
        cutout_sci, cutout_line, cutout_wht, cutout_dq,
        cutout_wht_drz, pixmap, pad_size, pad_half
    or None if the cutout falls outside the image boundary.
    """
    pad_half = int(np.ceil(cutout_size * np.sqrt(2) / 2)) + 1
    pad_size = 2 * pad_half + 1
    # cut window around the targeted emission line
    ny, nx = image.shape
    xs_int = int(np.round(xs))
    ys_int = int(np.round(ys))
    x_lo = xs_int - pad_half
    y_lo = ys_int - pad_half
    x_hi = xs_int + pad_half + 1
    y_hi = ys_int + pad_half + 1

    if x_lo < 0 or x_hi > nx or y_lo < 0 or y_hi > ny:
        return None

    cutout_sci  = image       [y_lo:y_hi, x_lo:x_hi].copy()
    cutout_line = emline      [y_lo:y_hi, x_lo:x_hi].copy()
    cutout_wht  = weight      [y_lo:y_hi, x_lo:x_hi].copy()
    cutout_dq   = data_quality[y_lo:y_hi, x_lo:x_hi].copy()

    cutout_wht_drz = np.nan_to_num(cutout_wht, nan=0.0, posinf=0.0, neginf=0.0)
    cutout_wht_drz[cutout_dq % 2 == 1] = 0.0

    # Build pixel map: each padded cutout pixel (ix, iy) maps to sky via the
    # *undispersed* position (x0 + offset, y0 + offset), eliminating sub-pixel
    # rounding error in xs_int = round(x0 + dx_line).
    iy_arr, ix_arr = np.mgrid[0:pad_size, 0:pad_size]
    gx_undis = (x0 + (ix_arr - pad_half)).ravel()
    gy_undis = (y0 + (iy_arr - pad_half)).ravel()
    # map the cutout pixel grid to sky coordinates at direct image position (x0, y0)
    sky_coords = wcs_di.all_pix2world(np.column_stack([gx_undis, gy_undis]), 0)
    # map the sky coordinates at direct image to output pixel coordinates
    out_pix    = out_wcs.all_world2pix(sky_coords, 0)
    pixmap = np.dstack([
        out_pix[:, 0].reshape(pad_size, pad_size),
        out_pix[:, 1].reshape(pad_size, pad_size),
    ])

    return dict(
        cutout_sci=cutout_sci, cutout_line=cutout_line,
        cutout_wht=cutout_wht, cutout_dq=cutout_dq,
        cutout_wht_drz=cutout_wht_drz, pixmap=pixmap,
        pad_size=pad_size, pad_half=pad_half,
    )


def _process_emline_frame(grism_fn, POM_fn, source_id, wave_line_obs,
                           grism_conf, filter, out_wcs, cutout_size,
                           psf_oversample, drz_psf, psf_meta, pixfrac, pixscale_ratio):
    """
    Drizzle helper for :func:`extract_2d_emline_worker_drizzle`.

    For a single grism frame, this function:

    1. Reads the POM catalog to get the undispersed source position ``(x0, y0)``.
    2. Loads the grism SCI, emline, DQ, and weight arrays.
    3. Calls ``get_position_ang_dispang_at_wave`` to obtain the predicted emission-line
       pixel position ``(xs, ys)`` and local dispersion angle ``theta`` in the
       sky-aligned output frame.
    4. Calls :func:`_build_cutout_and_pixmap` to produce a padded cutout and pixel map.
    5. Drizzles the SCI and LINE cutouts separately onto a ``cutout_size × cutout_size``
       sky-aligned output grid.
    6. Calls :func:`_compute_grism_psf_frame` and drizzles the resulting PSF into the
       shared ``drz_psf[mp_key]`` accumulator.

    Mutates ``drz_psf`` and ``psf_meta`` in place (PSF accumulation across frames).
    Returns ``None`` if the frame should be skipped (filter mismatch, source absent
    in POM catalog, or predicted position out of bounds).

    Returns
    -------
    result : dict or None
        Keys: ``cutout, xs, ys, theta, module, pupil, path, effexptm, gs_v3pa,
        mp_key, frame_sci, frame_line, frame_wht``.
        ``cutout`` is a tuple ``(sci, line, wht, dq)`` of 2-D arrays from the
        un-drizzled padded cutout.  ``frame_sci``, ``frame_line``, and ``frame_wht``
        are the drizzled outputs on the sky-aligned grid.
    """
    POM_cat    = ascii.read(POM_fn)
    primary_hd = fits.getheader(grism_fn)
    sci_hd     = fits.getheader(grism_fn, "sci")
    wcs_grism = WCS(sci_hd)

    _filter = primary_hd["FILTER"]
    module  = primary_hd["MODULE"]
    pupil   = primary_hd["PUPIL"][-1]

    if _filter != filter:
        return None

    image = fits.getdata(grism_fn, "sci")
    try:
        emline = fits.getdata(grism_fn, "emline")
    except KeyError:
        emline = image
    data_quality = fits.getdata(grism_fn, "dq")

    weight_path = grism_fn.replace("lv1.5.fits", "wht.fits")
    if os.path.isfile(weight_path):
        weight = fits.getdata(weight_path)
    else:
        weight = fits.getdata(grism_fn, "err")
        weight[weight == 0] = np.nan
        weight = weight ** -2

    item_POM = POM_cat[POM_cat["Index"] == source_id]
    if len(item_POM) == 0:
        return None
    x0 = float(item_POM["pixel_x"][0])
    y0 = float(item_POM["pixel_y"][0])

    ### Load dispersion solutions
    xs, ys, theta = get_position_ang_dispang_at_wave(
        x0, y0, wave_line_obs, grism_conf, module, pupil, 
        wcs_transform=[wcs_grism, out_wcs], velosys=sci_hd["velosys"])
    if (xs < VALID_X_LEFT) or (xs > VALID_X_RIGHT) or (ys < VALID_Y_BOTTOM) or (ys > VALID_Y_TOP):
        print(" >> [emline] ID%s %s%s predicted position (%.1f,%.1f) out of bounds; skip" %
              (source_id, module, pupil, xs, ys))
        return None

    ### Take cutout of emission line at the predicted position in grism image.
    cutout_result = _build_cutout_and_pixmap(
        image, emline, weight, data_quality,
        x0, y0, xs, ys, cutout_size * 1.5, wcs_grism, out_wcs,
    )
    if cutout_result is None:
        print(" >> [emline] ID%s %s%s cutout out of bounds at (%.1f,%.1f); skip" %
              (source_id, module, pupil, xs, ys))
        return None

    cutout_sci     = cutout_result['cutout_sci']
    cutout_line_   = cutout_result['cutout_line']
    cutout_wht     = cutout_result['cutout_wht']
    cutout_dq      = cutout_result['cutout_dq']
    cutout_wht_drz = cutout_result['cutout_wht_drz']
    pixmap         = cutout_result['pixmap']
    pad_size       = cutout_result['pad_size']

    mp_key   = "%s%s" % (module, pupil)
    effexptm = max(float(primary_hd.get("EFFEXPTM", 1.0)), 1e-6)
    gs_v3pa  = float(np.deg2rad(primary_hd.get("GS_V3_PA", np.nan)))

    if mp_key not in drz_psf:
        drz_psf[mp_key] = Drizzle(kernel="square", 
            out_shape=(cutout_size*psf_oversample, cutout_size*psf_oversample))
    out_wcs_psf = out_wcs.deepcopy()
    out_wcs_psf.wcs.cd = out_wcs_psf.wcs.cd / psf_oversample
    out_wcs_psf.array_shape = [out_wcs_psf.array_shape[0] * psf_oversample,
                                out_wcs_psf.array_shape[1] * psf_oversample]
    out_wcs_psf.wcs.set()

    drz_f_sci  = Drizzle(kernel="square", out_shape=(cutout_size, cutout_size))
    drz_f_line = Drizzle(kernel="square", out_shape=(cutout_size, cutout_size))
    drz_f_sci.add_image(
        cutout_sci.astype(np.float64), exptime=effexptm,
        pixmap=pixmap, weight_map=cutout_wht_drz.astype(np.float64),
        pixfrac=pixfrac, in_units="cps", scale=pixscale_ratio
    )
    drz_f_line.add_image(
        cutout_line_.astype(np.float64), exptime=effexptm,
        pixmap=pixmap, weight_map=cutout_wht_drz.astype(np.float64),
        pixfrac=pixfrac, in_units="cps", scale=pixscale_ratio
    )
    frame_sci  = drz_f_sci.out_img.copy()
    frame_line = drz_f_line.out_img.copy()
    frame_wht  = drz_f_sci.out_wht.copy()
    frame_sci [frame_wht == 0] = np.nan
    frame_line[frame_wht == 0] = np.nan

    # PSF model for this frame — drizzled into the shared accumulator
    psf_frame, psf_over_used, psf_ver = _compute_grism_psf_frame(
        primary_hd=primary_hd, xs=xs, ys=ys,
        filter_name=filter, oversample=psf_oversample, fov_pixels=pad_size,
    )
    if mp_key not in psf_meta:
        psf_meta[mp_key] = (psf_ver, psf_over_used)

    N_psf_over  = psf_frame.shape[0]
    iy_p, ix_p  = np.mgrid[0:N_psf_over, 0:N_psf_over]
    center_over = (N_psf_over - 1) / 2.0
    gx_psf = (x0 + (ix_p - center_over) / psf_over_used).ravel()
    gy_psf = (y0 + (iy_p - center_over) / psf_over_used).ravel()

    sky_psf  = wcs_grism.all_pix2world(np.column_stack([gx_psf, gy_psf]), 0)
    opix_psf = out_wcs_psf.all_world2pix(sky_psf, 0)
    pixmap_psf = np.dstack([
        opix_psf[:, 0].reshape(N_psf_over, N_psf_over),
        opix_psf[:, 1].reshape(N_psf_over, N_psf_over),
    ])
    drz_psf[mp_key].add_image(
        psf_frame.astype(np.float64), exptime=effexptm,
        pixmap=pixmap_psf, weight_map=np.ones_like(psf_frame, dtype=np.float64),
        pixfrac=pixfrac, in_units="cps", scale=pixscale_ratio
    )

    return dict(
        cutout=(cutout_sci, cutout_line_, cutout_wht, cutout_dq),
        xs=xs, ys=ys, theta=theta,
        module=module, pupil=pupil, path=grism_fn,
        effexptm=float(primary_hd.get("EFFEXPTM", np.nan)),
        gs_v3pa=gs_v3pa, mp_key=mp_key,
        frame_sci=frame_sci, frame_line=frame_line, frame_wht=frame_wht,
    )


def _sigma_clip_weighted_coadd(frame_list_sci, frame_list_line, frame_list_wht):
    """
    Sigma-clip (σ=2.5, 5 iters) across frames then compute inverse-variance
    weighted mean.  Returns (coadd_sci, coadd_line, coadd_wht, coadd_cov).
    """
    arr_sci  = np.array(frame_list_sci,  dtype=np.float64)
    arr_line = np.array(frame_list_line, dtype=np.float64)
    arr_wht  = np.array(frame_list_wht,  dtype=np.float64)

    if arr_sci.shape[0] > 2:
        sc = SigmaClip(sigma=2.5, maxiters=5)
        arr_sci  = sc(arr_sci,  axis=0, masked=False)
        arr_line = sc(arr_line, axis=0, masked=False)

    arr_wht[np.isnan(arr_sci)]  = 0.0
    arr_wht[np.isnan(arr_line)] = 0.0

    wht_sum    = np.nansum(arr_wht, axis=0)
    coadd_sci  = np.where(wht_sum > 0,
                          np.nansum(arr_sci  * arr_wht, axis=0) / wht_sum, np.nan)
    coadd_line = np.where(wht_sum > 0,
                          np.nansum(arr_line * arr_wht, axis=0) / wht_sum, np.nan)
    coadd_wht  = wht_sum
    coadd_cov  = np.int8(np.sum(arr_wht > 0, axis=0))
    return coadd_sci, coadd_line, coadd_wht, coadd_cov


def _build_emline_coadd_hdul(
    source_id, source_ra, source_dec,
    filter_name, module_mp, pupil_mp,
    name_line, wave_line_obs, eml_snr, z_grism,
    pixscale_out,
    coadd_sci, coadd_line, coadd_wht, coadd_cov,
    coadd_psf,
    n_coadd, mean_theta, diff_theta, mean_effexptm, mean_gs_v3pa,
    frame_pa=0.0,
    out_wcs_header=None,
    psf_meta_mp=None,
    wave_1=None,
    d_wave=None,
    cutout_size=None,
    chiral=None,
):
    """Build the output HDUList for one module-pupil emission-line coadd.

    All angle parameters (mean_theta, diff_theta, mean_gs_v3pa, frame_pa) must
    be in radians; they are stored as-is in FITS headers with unit comment "(rad)".

    Returns an ``HDUList`` with extensions:
    ``PRIMARY``, ``SPEC2D``, ``WHT2D``, ``COV2D``, ``LINE2D``, ``PSF2D``.

    Parameters
    ----------
    source_id : int or str
        Source identifier written to the ``ID`` header keyword.
    source_ra, source_dec : float
        Source sky position in degrees.
    filter_name : str
        Grism filter name (e.g. ``'F444W'``).
    module_mp, pupil_mp : str
        Detector module (``'A'`` or ``'B'``) and pupil (``'R'`` or ``'C'``).
    name_line : str
        Name of the targeted emission line (e.g. ``'Ha'``).
    wave_line_obs : float
        Observed wavelength of the emission line in µm.
    eml_snr : float
        Best-fit S/N of the line across all frames.
    z_grism : float
        Spectroscopic redshift from the grism.
    pixscale_out : float
        Output pixel scale in arcsec/pixel.
    coadd_sci, coadd_line : 2-D ndarray
        Co-added SCI and continuum-subtracted LINE images.
    coadd_wht : 2-D ndarray
        Co-added inverse-variance weight map.
    coadd_cov : 2-D ndarray
        Integer coverage map (number of frames covering each pixel).
    coadd_psf : 2-D ndarray
        Co-added (and normalised) PSF model image.
    n_coadd : int
        Number of frames included in the coadd.
    mean_theta : float
        Mean dispersion angle w.r.t. the cutout +X axis (rad, CCW from +X).
    diff_theta : float
        Frame-to-frame variation of the dispersion angle (rad).
    mean_effexptm : float
        Total effective exposure time (s).
    mean_gs_v3pa : float
        Mean V3 position angle across frames (rad).
    frame_pa : float, optional
        Position angle of the cutout +X axis (northward from west, rad).
        Default is 0.0 (sky-aligned drizzle case).
    out_wcs_header : dict or None
        Sky-aligned WCS key-value pairs (drizzle case); ``None`` for simple coadd.
    psf_meta_mp : (str, int) or None
        ``(version_string, oversample)`` tuple used to populate PSF extension
        headers.
    wave_1 : float or None
        Wavelength of the first cutout pixel in µm (simple case only).
    d_wave : float or None
        Wavelength step in µm/pix (simple case only).
    cutout_size : int or None
        Side length of the cutout in pixels; stored as ``APERTURE`` (simple case only).
    chiral : int or None
        Chirality convention: +1 right-handed, -1 left-handed (simple case only).
    """
    # --- Primary HDU ---------------------------------------------------------
    hdu = fits.PrimaryHDU()
    hdu.header["ID"]       = (source_id,     "Source ID")
    hdu.header["RA0"]      = (source_ra,     "Source RA (deg)")
    hdu.header["DEC0"]     = (source_dec,    "Source DEC (deg)")
    hdu.header["FILTER"]   = (filter_name,   "Filter name")
    hdu.header["MODULE"]   = (module_mp,     "Detector module (A or B)")
    hdu.header["PUPIL"]    = (pupil_mp,      "Pupil (R=GRISMR, C=GRISMC)")
    hdu.header["LINENAME"] = (name_line,     "Target emission line")
    hdu.header["LINEWAVE"] = (wave_line_obs, "Observed wavelength of line (um)")
    hdu.header["LINESNR"]  = (eml_snr,       "Bestfit S/N of the line from all frames")
    hdu.header["ZGRISM"]   = (z_grism,       "Spectroscopic redshift from grism")
    hdu.header["N_COADD"]  = (n_coadd,       "Number of coadded frames")
    hdu.header["DISPANG"]  = (mean_theta,
                              "Mean dispersion angle w.r.t. the cutout X axis (rad, CCW from +X)")
    hdu.header["DISPANGW"] = (diff_theta,
                              "Dispersion angle variation among frames (rad)")
    hdu.header["FRAME_PA"] = (frame_pa,
                              "Position angle of the cutout X axis (northward from west, rad)")
    hdu.header["GS_V3_PA"] = (mean_gs_v3pa,
                              "Mean V3 position angle (rad)")
    hdu.header["EFFEXPTM"] = (mean_effexptm, "Total effective exposure time (s)")
    hdu.header["PIXSCALE"] = (pixscale_out,  "Spatial pixel scale (arcsec/pix)")
    if wave_1 is not None:
        hdu.header["WAVE_1"] = (wave_1, "Wavelength of first cutout pixel (um)")
    if d_wave is not None:
        hdu.header["D_WAVE"] = (d_wave, "Wavelength step (um/pix)")
    hdu.header["AUTHOR"]   = ("Jiachuan Xu", "Author")
    hdu.header["TIME"]     = (
        time.strftime("%Y/%m/%d %H:%M:%S", time.localtime()), "Creation time"
    )
    if out_wcs_header is not None:
        for key, val in out_wcs_header.items():
            hdu.header[key] = val

    # --- Image extensions ----------------------------------------------------
    frame_comment = ("+X = west (+RA), +Y = north (+DEC)" if out_wcs_header is not None
                     else "+X = dispersion (+wavelength), +Y = cross-dispersion")

    hdu_sci  = fits.ImageHDU(coadd_sci.astype(np.float32),  name="SPEC2D")
    hdu_wht  = fits.ImageHDU(coadd_wht.astype(np.float32),  name="WHT2D")
    hdu_cov  = fits.ImageHDU(coadd_cov.astype(np.int8),     name="COV2D")
    hdu_line = fits.ImageHDU(coadd_line.astype(np.float32), name="LINE2D")

    for _hdu in [hdu_sci, hdu_line]:
        _hdu.header["BUNIT"]    = ("DN/s", "Brightness unit (no sensitivity correction)")
        _hdu.header["COMMENT"]  = frame_comment
        _hdu.header["PIXSCALE"] = (pixscale_out, "Spatial pixel scale (arcsec/pix)")
        _hdu.header["DISPANG"]  = (float(mean_theta),
                                   "Mean dispersion angle w.r.t. the cutout X axis (rad, CCW from +X)")
        _hdu.header["DISPANGW"] = (float(diff_theta),
                                   "Dispersion angle variation among frames (rad)")
        _hdu.header["FRAME_PA"] = (float(frame_pa),
                                   "Position angle of the cutout X axis (northward from west, rad)")
        _hdu.header["EMLFITSN"] = (eml_snr, "Bestfit S/N of the line from all frames")
        _hdu.header["EMISSION"] = (name_line, "Target emission line")
        if cutout_size is not None:
            _hdu.header["APERTURE"] = (cutout_size, "Cutout size (pixels)")
        if wave_1 is not None:
            _hdu.header["WAVE_1"] = (wave_1, "Wavelength of first cutout pixel (um)")
        if d_wave is not None:
            _hdu.header["D_WAVE"] = (d_wave, "Wavelength step (um/pix)")
        if chiral is not None:
            _hdu.header["CHIRAL"] = (chiral, "Chirality: +1 right-handed, -1 left-handed")

    # --- PSF extension -------------------------------------------------------
    _psf_ver, _psf_over = psf_meta_mp if psf_meta_mp is not None else ("unknown", 1)
    hdu_psf = fits.ImageHDU(coadd_psf.astype(np.float32), name="PSF2D")
    hdu_psf.header["BUNIT"]    = ("",           "Normalised PSF (sums to 1)")
    hdu_psf.header["COMMENT"]  = frame_comment
    hdu_psf.header["NAXIS"]  = 2
    hdu_psf.header["NAXIS1"] = int(coadd_psf.shape[1])
    hdu_psf.header["NAXIS2"] = int(coadd_psf.shape[0])
    hdu_psf.header["MODULE"] = (module_mp,     "Detector module (A or B)")
    hdu_psf.header["WAVELEN"] = (wave_line_obs*1e-6, "Wavelength at which PSF is computed (m)")
    hdu_psf.header["CHANNEL"] = ("Long",   "Detector channel (Short or Long)")
    hdu_psf.header["PILIN"] = False
    hdu_psf.header["FILTER"] = (filter_name,  "Filter used for PSF calculation")
    hdu_psf.header["PIXELSCL"] = (pixscale_out, "Output pixel scale (arcsec/pix)")
    hdu_psf.header["OVERSAMP"] = (_psf_over,    "PSF super-sampling factor (native pix)")
    hdu_psf.header["DET_SAMP"] = (_psf_over,    "PSF super-sampling factor (native pix)")
    hdu_psf.header["DET_NAME"] = ("NRC%s5" % module_mp, "Detector name passed to stpsf")
    hdu_psf.header["DIFFLMT"] = (1.22 * wave_line_obs*1e-6/6.5 * 206265,
                                    "Diffraction limit FWHM (arcsec) at the line wavelength")
    hdu_psf.header["FOV"] = (pixscale_out * coadd_psf.shape[0], "PSF field of view (arcsec)")
    hdu_psf.header["N_COADD"] = (n_coadd,      "Number of PSF frames coadded")
    if out_wcs_header is not None:
        for key, val in out_wcs_header.items():
            hdu_psf.header[key] = val

    return fits.HDUList([hdu, hdu_sci, hdu_wht, hdu_cov, hdu_line, hdu_psf])


def extract_2d_emline_worker_drizzle(
    grism_idx_per_source: list,
    POM_path_per_source: list,
    all_v1p5_list: list,
    source_item: "Table.Row",
    grism_conf: "GrismConf",
    filter: str,
    extraction_dir: str,
    cutout_size: int = 51,
    cutout_pixscale: float = 0.031,
    pixfrac: float = 1.0,
    psf_oversample: int = 4,
) -> None:
    """
    Extract drizzle-coadded 2-D emission-line cutouts for a single source.

    For each grism frame where the source yields a spectrum, this function:

    1. Predicts the pixel position of the target emission line from the
       source redshift and the grism dispersion model.
    2. Takes an ``cutout_size × cutout_size`` cutout (in DN/s, no sensitivity
       correction) from both the continuum-included (SCI) and continuum-subtracted
       (LINE) grism images.
    3. Computes the local dispersion angle θ = angle(dX/dλ + j·dY/dλ) in the
       sky-aligned output frame (+X west, +Y north).
    4. Maps each cutout pixel to the sky-aligned output frame by first removing
       the spectral offset (undispersed position), applying the grism WCS, and
       then reprojecting via the output WCS.
    5. Drizzles each per-frame cutout individually onto the output grid, then
       sigma-clips (σ=2.5, 5 iterations) across frames at each output pixel
       and computes an inverse-variance weighted mean, separately for each
       module-pupil combination.  Output pixel scale =
       ``cutout_pixscale`` arcsec/pix; output grid size =
       ``cutout_size`` pixels
       on a side; drop-size fraction = ``pixfrac``.
    6. Renders a super-sampled PSF model for each frame via
       :func:`_compute_grism_psf_frame` (using ``stpsf`` with in-flight OPD
       when available, falling back to a Gaussian approximation) and drizzles
       it into the same output grid.  The coadded PSF model is saved as the
       ``PSF2D`` extension in the per-module-pupil coadd FITS file.

    Results are saved as:

    * ``emline_2d_{filter}_ID{id}_{line}_all.fits``   – all per-frame cutouts.
    * ``emline_2d_{filter}_ID{id}_{line}_{MP}coadd_drz.fits`` – drizzled coadd
      per module-pupil (extensions: SPEC2D, WHT2D, COV2D, LINE2D, PSF2D).

    Parameters
    ----------
    grism_idx_per_source:
        Indices into ``all_v1p5_list`` for frames where the source is observable.
    POM_path_per_source:
        Paths to the POM-applied source catalogs for those frames.
    all_v1p5_list:
        List of all level-1.5 grism FITS file paths.
    source_item:
        Single row from the source catalog; must contain ``ID``, ``RA``, ``DEC``,
        ``z_grism``, ``fit_line_SN``, and ``name_line_exp`` columns.
    grism_conf:
        :class:`GrismConf` instance with dispersion and sensitivity calibration.
    filter:
        Grism filter name (e.g. ``'F356W'``).
    extraction_dir:
        Output directory for the drizzled FITS files.
    cutout_size:
        Side length of the drizzled output image in pixels.  Each per-frame
        cutout is padded before drizzling; the final output grid is exactly
        ``cutout_size × cutout_size`` pixels at ``cutout_pixscale`` arcsec/pixel.
    cutout_pixscale:
        Pixel scale of the drizzled output image in arcsec/pixel.
    pixfrac:
        Drizzle drop-shrunken scale factor (0 < pixfrac ≤ 1).
    psf_oversample:
        Super-sampling factor for the stpsf PSF model.  The PSF is rendered at
        ``psf_oversample × native_pixel_scale`` before being drizzled onto the
        output grid.  Higher values give a more accurate PSF at the cost of
        longer stpsf computation time.  Ignored in the Gaussian fallback, which
        always produces a native-scale PSF (oversample_used = 1).
    """
    if len(grism_idx_per_source) == 0:
        return

    source_id    = source_item["ID"]
    source_ra    = float(source_item["RA"])
    source_dec   = float(source_item["DEC"])
    eml_snr      = float(source_item["fit_line_SN"])
    source_coord = SkyCoord(source_ra, source_dec, unit=(u.deg, u.deg))

    try:
        name_line = str(source_item["name_line_exp"]).strip()
        z_grism   = float(source_item["z_grism"])
    except (KeyError, ValueError, TypeError) as exc:
        print(" >> [emline] source %s: missing z_grism / name_line_exp (%s); skip" %
              (source_id, exc))
        return

    if name_line not in EML_LAB:
        print(" >> [emline] source %s: unknown line '%s'; skip" % (source_id, name_line))
        return

    wave_line_rest = EML_LAB[name_line]
    wave_line_obs  = wave_line_rest * (1.0 + z_grism)

    if not (grism_conf.WRANGE[0] <= wave_line_obs <= grism_conf.WRANGE[1]):
        print(" >> [emline] source %s: line %.4f µm outside filter range; skip" %
              (source_id, wave_line_obs))
        return

    line_tag = "".join(c if c.isalnum() else "_" for c in name_line)

    # Output WCS: sky-aligned, +X=west, +Y=north.
    out_wcs, out_wcs_header = _make_sky_wcs(source_ra, source_dec, cutout_pixscale, cutout_size)
    pixscale_ratio = cutout_pixscale / NIRCAM_LW_PIXSCALE

    # Per-frame accumulators
    cutout_list   = []
    xs_list       = []
    ys_list       = []
    theta_list    = []
    module_list   = []
    pupil_list    = []
    paths_list    = []
    effexptm_list = []
    gs_v3pa_list  = []

    frames_sci:  dict[str, list] = {}
    frames_line: dict[str, list] = {}
    frames_wht:  dict[str, list] = {}
    drz_psf:     dict[str, Drizzle] = {}
    psf_meta:    dict[str, tuple]   = {}

    for j, POM_fn in enumerate(POM_path_per_source):
        grism_fn = all_v1p5_list[grism_idx_per_source[j]]
        result = _process_emline_frame(
            grism_fn, POM_fn, source_id, wave_line_obs,
            grism_conf, filter, out_wcs, cutout_size,
            psf_oversample, drz_psf, psf_meta, pixfrac, pixscale_ratio,
        )
        if result is None:
            continue

        cutout_list.append(result['cutout'])
        xs_list.append(result['xs'])
        ys_list.append(result['ys'])
        theta_list.append(result['theta'])
        module_list.append(result['module'])
        pupil_list.append(result['pupil'])
        paths_list.append(result['path'])
        effexptm_list.append(result['effexptm'])
        gs_v3pa_list.append(result['gs_v3pa'])

        mp_key = result['mp_key']
        if mp_key not in frames_sci:
            frames_sci[mp_key]  = []
            frames_line[mp_key] = []
            frames_wht[mp_key]  = []
        frames_sci[mp_key].append(result['frame_sci'])
        frames_line[mp_key].append(result['frame_line'])
        frames_wht[mp_key].append(result['frame_wht'])

    if len(cutout_list) == 0:
        print(" >> [emline] ID%s: no valid cutouts found" % source_id)
        return

    all_fn = os.path.join(
        extraction_dir,
        "emline_2d_%s_ID%s_%s_all.fits" % (filter, source_id, line_tag),
    )
    store_all_2d_emline(
        cutout_list, xs_list, ys_list, theta_list,
        module_list, pupil_list, paths_list, effexptm_list, gs_v3pa_list,
        output=all_fn, coord=source_coord, source_item=source_item,
        grism_filter=filter, wave_line=wave_line_obs, name_line=name_line,
    )
    print(" >> [emline] ID%s saved per-frame: %s" % (source_id, all_fn))

    for mp_key in frames_sci:
        module_mp = mp_key[0]
        pupil_mp  = mp_key[1]

        idx_mp = [i for i, (m, p) in enumerate(zip(module_list, pupil_list))
                  if m == module_mp and p == pupil_mp]
        if not idx_mp:
            continue

        mean_theta    = float(np.nanmean([theta_list[i]    for i in idx_mp]))
        diff_theta    = float(np.nanmax([theta_list[i]     for i in idx_mp]) -
                              np.nanmin([theta_list[i]     for i in idx_mp]))
        mean_effexptm = float(np.nansum([effexptm_list[i]  for i in idx_mp]))
        mean_gs_v3pa  = float(np.nanmean([gs_v3pa_list[i]  for i in idx_mp]))

        coadd_sci, coadd_line, coadd_wht, coadd_cov = _sigma_clip_weighted_coadd(
            frames_sci[mp_key], frames_line[mp_key], frames_wht[mp_key],
        )

        coadd_psf_raw = drz_psf[mp_key].out_img
        psf_total     = coadd_psf_raw.sum()
        coadd_psf     = coadd_psf_raw / psf_total if psf_total > 0 else coadd_psf_raw

        hdul_out = _build_emline_coadd_hdul(
            source_id=source_id, source_ra=source_ra, source_dec=source_dec,
            filter_name=filter, module_mp=module_mp, pupil_mp=pupil_mp,
            name_line=name_line, wave_line_obs=wave_line_obs,
            eml_snr=eml_snr, z_grism=z_grism,
            pixscale_out=cutout_pixscale,
            coadd_sci=coadd_sci, coadd_line=coadd_line,
            coadd_wht=coadd_wht, coadd_cov=coadd_cov,
            coadd_psf=coadd_psf,
            n_coadd=len(idx_mp),
            mean_theta=mean_theta, diff_theta=diff_theta,
            mean_effexptm=mean_effexptm,
            mean_gs_v3pa=mean_gs_v3pa,
            frame_pa=0.0,
            out_wcs_header=out_wcs_header,
            psf_meta_mp=psf_meta.get(mp_key, ("unknown", psf_oversample)),
        )
        coadd_fn = os.path.join(
            extraction_dir,
            "emline_2d_%s_ID%s_%s_%scoadd_drz.fits" % (filter, source_id, line_tag, mp_key),
        )
        hdul_out.writeto(coadd_fn, overwrite=True)
        print(" >> [emline] ID%s %s coadd → %s" % (source_id, mp_key, coadd_fn))


def extract_2d_emline_worker_simple(
    source_item: "Table.Row",
    grism_conf: "GrismConf",
    filter: str,
    extraction_dir: str,
    cutout_size: int = 51,
    psf_oversample: int = 4,
) -> None:
    """
    Extract simple-coadded 2-D emission-line cutouts for a single source.

    This function directly take cutouts from the coadded 2D spectra in
    :func:`extract_2d_spec_worker` method, without performing drizzle re-sampling,
    and assuming a unified charility of the coadded 2D spectra of +1. Therefore, 
    when measuring kinematic lensing, the galaxy model will be rotated into this 
    spectra frame to fit the spectrum data. This means when running with cutouts
    produced by this function, KL always assume the dispersion angle is almost 0 
    (disperse along +X direction) and the spectra are right-handed. 
    For each source, this function will:

    1. Find out the relavant coadded 2D spectra from `extraction_dir` by matching
       the filename convention ``spec_2d_{FILTER}_ID{OBJID}_{MODULE}{PUPIL}coadd.fits`` 
    2. Given the input ``source_item``, identify the target emission line and its 
       observed wavelength, and check if the line falls within the wavelength range
       of the coadded 2D spectra. If so, take a ``cutout_size × cutout_size`` cutout
       from both the continuum-included (SCI) and continuum-subtracted (LINE) extensions
       of the coadded 2D spectra, centred at the predicted pixel position of the target 
       emission line.
       Note that the cutouts are taken in the DETECTOR frame, NOT equatorial coord frame. 
       Therefore the +X axis has a non-trivial position angle in the sky.
    3. Computes the local dispersion angle θ = angle(dX/dλ + j·dY/dλ) in the DETECTOR 
       frame (``DISPANG``), and its standard deviation. Also compute the position angle 
       (``FRAME_PA``) of the +X axis in EQUATORIAL frame (i.e. north up, east left).
    4. Renders a super-sampled PSF model for each frame via
       :func:`_compute_grism_psf_frame` (using ``stpsf`` with in-flight OPD
       when available, falling back to a Gaussian approximation) and combines
       it into a super-sampled output grid.  The coadded PSF model is saved as the
       ``PSF2D`` extension in the per-module-pupil coadd FITS file.

    Results are saved as:

    * ``emline_2d_{filter}_ID{id}_{line}_{MP}coadd_simple.fits`` – simple coadd
      per module-pupil (extensions: SPEC2D, WHT2D, COV2D, LINE2D, PSF2D).

    Parameters
    ----------
    source_item:
        Single row from the source catalog; must contain ``ID``, ``RA``, ``DEC``,
        ``z_grism``, ``fit_line_SN``, and ``name_line_exp`` columns.
    grism_conf : GrismConf
        Grism configuration object containing dispersion calibration, used to
        compute the predicted emission-line position in each coadded 2D spectrum.
    filter:
        Grism filter name (e.g. ``'F356W'``).
    extraction_dir:
        Output directory for the simple-coadd FITS files.
    cutout_size:
        Side length of the output cutout in native detector pixels (must be odd
        for a symmetric aperture).  The output is exactly ``cutout_size ×
        cutout_size`` pixels at the native NIRCam LW pixel scale — no drizzle
        resampling is performed.
    psf_oversample:
        Super-sampling factor for the stpsf PSF model.  The PSF is rendered at
        ``psf_oversample × native_pixel_scale`` before being drizzled onto the
        output grid.  Higher values give a more accurate PSF at the cost of
        longer stpsf computation time.  Ignored in the Gaussian fallback, which
        always produces a native-scale PSF.
    """
    # --- Source info ----------------------------------------------------------
    source_id  = source_item["ID"]
    source_ra  = float(source_item["RA"])
    source_dec = float(source_item["DEC"])
    eml_snr    = float(source_item["fit_line_SN"])

    # --- Emission line info ---------------------------------------------------
    try:
        name_line = str(source_item["name_line_exp"]).strip()
        z_grism   = float(source_item["z_grism"])
    except (KeyError, ValueError, TypeError) as exc:
        print(" >> [emline_simple] source %s: missing z_grism / name_line_exp (%s); skip" %
              (source_id, exc))
        return

    if name_line not in EML_LAB:
        print(" >> [emline_simple] source %s: unknown line '%s'; skip" % (source_id, name_line))
        return

    wave_line_rest = EML_LAB[name_line]
    wave_line_obs  = wave_line_rest * (1.0 + z_grism)

    line_tag = "".join(c if c.isalnum() else "_" for c in name_line)

    # --- Loop over module-pupil combinations ---------------------------------
    for module in ["A", "B"]:
        for pupil in ["R", "C"]:
            mp_key   = "%s%s" % (module, pupil)
            coadd_fn = os.path.join(
                extraction_dir,
                "spec_2d_%s_ID%s_%scoadd.fits" % (filter, source_id, mp_key),
            )
            if not os.path.isfile(coadd_fn):
                continue

            with fits.open(coadd_fn) as hdul:
                sci_hd    = hdul["SPEC2D"].header
                sci_data  = hdul["SPEC2D"].data.copy()    # (N_spat, N_wave)
                wht_data  = hdul["WHT2D"].data.copy()
                cov_data  = hdul["COV2D"].data.copy()
                line_data = hdul["LINE2D"].data.copy()
                stats_tb  = Table(hdul["STATS"].data)

            # --- Wavelength calibration and range check -----------------------
            wave_1  = sci_hd["wave_1"]
            d_wave  = sci_hd["d_wave"]
            n_spat, n_wave = sci_data.shape

            if wave_line_obs < wave_1 or wave_line_obs > wave_1 + (n_wave - 1) * d_wave:
                print(" >> [emline_simple] ID%s %s: line %.4f µm outside coadd range "
                      "[%.4f, %.4f]; skip" %
                      (source_id, mp_key, wave_line_obs,
                       wave_1, wave_1 + (n_wave - 1) * d_wave))
                continue

            # --- Cutout centre ------------------------------------------------
            # Emission-line column (wavelength axis, X = NAXIS1 direction).
            col_line   = (wave_line_obs - wave_1) / d_wave
            # Source is placed at the spatial centre of the coadded 2D spec.
            row_center = (n_spat - 1) / 2.0

            half    = (cutout_size - 1) // 2
            col_int = int(np.round(col_line))
            row_int = int(np.round(row_center))
            col_lo  = col_int - half
            col_hi  = col_int + half + 1   # exclusive → cutout_size columns
            row_lo  = row_int - half
            row_hi  = row_int + half + 1   # exclusive → cutout_size rows

            if col_lo < 0 or col_hi > n_wave or row_lo < 0 or row_hi > n_spat:
                print(" >> [emline_simple] ID%s %s: cutout out of bounds at "
                      "(col=%.1f, row=%.1f); skip" %
                      (source_id, mp_key, col_line, row_center))
                continue

            cutout_sci  = sci_data [row_lo:row_hi, col_lo:col_hi].copy()
            cutout_wht  = wht_data [row_lo:row_hi, col_lo:col_hi].copy()
            cutout_cov  = cov_data [row_lo:row_hi, col_lo:col_hi].copy()
            cutout_line = line_data[row_lo:row_hi, col_lo:col_hi].copy()

            # --- Dispersion angle and position angle --------------------------
            # In the 2D-spec frame, wavelength is always along +X by construction,
            # so the dispersion angle is almost zero with no frame-to-frame
            # variation.
            dispang_list = []
            xs_list, ys_list, date_beg_list = [], [], []
            for iexp in range(len(stats_tb)):
                x0 = stats_tb["x0"][iexp]
                y0 = stats_tb["y0"][iexp]
                xs, ys, theta = get_position_ang_dispang_at_wave(
                    x0, y0, wave_line_obs, grism_conf, module, pupil)
                if (xs < VALID_X_LEFT) or (xs > VALID_X_RIGHT) or (ys < VALID_Y_BOTTOM) or (ys > VALID_Y_TOP):
                    continue
                dispang_list.append(theta)
                xs_list.append(xs)
                ys_list.append(ys)
                date_beg_list.append(stats_tb["DATE-BEG"][iexp])
            n_coadd = len(dispang_list)
            if n_coadd == 0:
                print(" >> [emline_simple] ID%s %s: no valid emission line cutout" %
                      (source_id, mp_key))
                continue
            theta_det     = float(np.nanmean(dispang_list))
            theta_det_std = float(np.nanstd(dispang_list))
            if not np.isfinite(theta_det):
                print(" >> [emline_simple] ID%s %s: invalid mean dispersion angle" %
                      (source_id, mp_key))
                print("    dispersion angles:", dispang_list)
                exit(1)

            # Position angle of the +X (dispersion direction) axis in equatorial frame 
            # (north up, east left), derived from GS_V3_PA:
            mean_gs_v3pa  = float(np.deg2rad(np.nanmean(stats_tb["GS_V3_PA"]%360.0)))
            mean_effexptm = float(np.nansum(stats_tb["EFFEXPTM"]))
            if mp_key == "AC":
                pa_x_eq = (mean_gs_v3pa + np.pi/2.) % (2*np.pi)
                theta_det = (theta_det - np.pi/2.) % (2*np.pi)
            elif mp_key == "BR":
                pa_x_eq = (mean_gs_v3pa + np.pi) % (2*np.pi)
                theta_det = (theta_det - np.pi) % (2*np.pi)
            elif mp_key == "BC":
                pa_x_eq = (mean_gs_v3pa + np.pi/2.) % (2*np.pi)
                theta_det = (theta_det - np.pi/2.) % (2*np.pi)
            else:   # "AR"
                pa_x_eq = mean_gs_v3pa

            # --- PSF model ----------------------------------------------------
            # Mean PSF from stpsf across frames, at the predicted line position.
            psf_over_raw_list = []
            for _xs, _ys, _d in zip(xs_list, ys_list, date_beg_list):
                psf_hd  = fits.Header()
                psf_hd["MODULE"]   = module
                psf_hd["DET_NAME"] = "NRC%s5" % module
                psf_hd["DATE-BEG"] = _d
                psf_over_raw, _, _ = _compute_grism_psf_frame(
                    primary_hd=psf_hd,
                    xs=_xs,
                    ys=_ys,
                    filter_name=filter,
                    oversample=psf_oversample,
                    fov_pixels=cutout_size,
                )
                psf_over_raw_list.append(psf_over_raw)
            psf_over_raw = np.array(psf_over_raw_list).mean(axis=0)
            psf_norm = psf_over_raw / psf_over_raw.sum()

            # --- Build and save output FITS file -------------------------------------
            hdul_out = _build_emline_coadd_hdul(
                source_id=source_id, source_ra=source_ra, source_dec=source_dec,
                filter_name=filter, module_mp=module, pupil_mp=pupil,
                name_line=name_line, wave_line_obs=wave_line_obs,
                eml_snr=eml_snr, z_grism=z_grism,
                pixscale_out=NIRCAM_LW_PIXSCALE,
                coadd_sci=cutout_sci, coadd_line=cutout_line,
                coadd_wht=cutout_wht, coadd_cov=cutout_cov,
                coadd_psf=psf_norm,
                n_coadd=n_coadd,
                mean_theta=theta_det,
                diff_theta=theta_det_std,
                mean_effexptm=mean_effexptm,
                mean_gs_v3pa=mean_gs_v3pa,
                frame_pa=pa_x_eq,
                wave_1=wave_1 + col_lo * d_wave,
                d_wave=d_wave,
                cutout_size=cutout_size,
                chiral=1,
            )
            out_fn = os.path.join(
                extraction_dir,
                "emline_2d_%s_ID%s_%s_%scoadd_simple.fits" % (
                    filter, source_id, line_tag, mp_key),
            )
            hdul_out.writeto(out_fn, overwrite=True)
            print(" >> [emline_simple] ID%s %s coadd → %s" % (source_id, mp_key, out_fn))

    return