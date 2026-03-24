"""
blending.py
-----------
Tools for detecting galaxy blending and estimating flux contamination in
JWST/HST imaging data.

Primary function
----------------
check_blending(ra, dec, images, ...) → dict

The function accepts one or more galaxy positions (RA, Dec) and one or more
image bands.  For each galaxy it:
  - cuts out a postage-stamp image around the target
  - estimates and subtracts the local background (sep)
  - detects and deblends all sources in the stamp (sep)
  - labels each pixel with its parent source via a segmentation map
  - measures the total flux inside a circular aperture at the target position
  - measures the contaminating flux (pixels NOT belonging to the target) inside
    the same aperture
  - returns:
      is_blended              – bool flag per galaxy (and per band if multiband)
      contamination_fraction  – (total_flux - target_flux) / total_flux
      nearest_sep_arcsec      – arcsec to the nearest other detected source
      n_blending_neighbors    – number of neighbors within blend_sep_arcsec

Dependencies
------------
  numpy, astropy, sep  (all available in the klpipe conda env)
  scipy  (optional; used only for a fallback segmentation method)

Example
-------
>>> import astropy.io.fits as fits
>>> from astropy.wcs import WCS
>>> from blending import check_blending
>>>
>>> hdul = fits.open("jades_f444w.fits")
>>> data = hdul[1].data
>>> wcs  = WCS(hdul[1].header)
>>>
>>> ra  = [53.1234, 53.1240]
>>> dec = [-27.5678, -27.5690]
>>>
>>> # single-band
>>> results = check_blending(ra, dec, (data, wcs))
>>>
>>> # multi-band
>>> images = {"F200W": (data_f200w, wcs_f200w),
...           "F444W": (data_f444w, wcs_f444w)}
>>> results = check_blending(ra, dec, images, multiband=True)
>>>
>>> print(results["is_blended"])          # shape (N,)
>>> print(results["contamination_fraction"])  # shape (N,) or (N, Nbands)
>>> print(results["nearest_sep_arcsec"])  # shape (N,)
"""

import sys
import warnings
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.nddata import Cutout2D, NoOverlapError
from astropy.wcs.utils import proj_plane_pixel_scales

import sep


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _sep_array(arr):
    """Return a C-contiguous, native-byte-order float64 array for sep."""
    arr = np.array(arr, dtype=np.float64)
    # sep requires native byte order
    if arr.dtype.byteorder not in ("=", "|", sys.byteorder):
        arr = arr.byteswap().newbyteorder()
    return np.ascontiguousarray(arr)


def _pixel_scale_arcsec(wcs_obj):
    """Return mean pixel scale in arcsec/pixel."""
    try:
        scales = proj_plane_pixel_scales(wcs_obj)        # degrees/pixel
        return float(np.mean(np.abs(scales))) * 3600.0
    except Exception:
        try:
            return abs(float(wcs_obj.wcs.cdelt[0])) * 3600.0
        except Exception:
            return 0.031  # JWST NIRCam SW default, arcsec/pixel


def _background_subtract(data, bkg_box_size, bkg_filtersize):
    """
    Estimate and subtract a 2-D background with sep.
    Returns (bkg_subtracted_image, background_rms_image).
    """
    box = max(4, bkg_box_size)
    box = min(box, min(data.shape) // 2)   # can't exceed half the stamp
    try:
        bkg = sep.Background(data, bw=box, bh=box,
                             fw=bkg_filtersize, fh=bkg_filtersize)
        return (data - bkg.back()).astype(np.float64), bkg.rms().astype(np.float64)
    except Exception:
        # Fallback: simple sigma-clipped stats
        from astropy.stats import sigma_clipped_stats
        _, med, std = sigma_clipped_stats(data, sigma=3.0)
        return (data - med).astype(np.float64), np.full_like(data, std, dtype=np.float64)


def _extract_sources(bkg_sub, rms, detection_sigma, npixels,
                     deblend_nthresh, deblend_cont):
    """
    Detect and deblend sources using sep.  Returns (objects, segmap).
    segmap is None if the installed sep version does not support it.
    Uses a scalar detection threshold (global RMS median × sigma) for
    compatibility across sep versions.
    """
    # sep.extract expects a scalar threshold in many builds; use the
    # global (median) RMS so that faint-source environments are handled
    # consistently across the stamp.
    thresh_scalar = float(detection_sigma * np.median(rms))

    try:
        objects, segmap = sep.extract(
            bkg_sub, thresh_scalar,
            minarea=npixels,
            deblend_nthresh=deblend_nthresh,
            deblend_cont=deblend_cont,
            segmentation_map=True,
        )
        return objects, segmap.astype(np.int32)
    except TypeError:
        # Either segmentation_map not supported, or thresh must be scalar
        try:
            objects = sep.extract(
                bkg_sub, thresh_scalar,
                minarea=npixels,
                deblend_nthresh=deblend_nthresh,
                deblend_cont=deblend_cont,
            )
            return objects, None
        except Exception:
            return np.array([]), None


def _find_target_label(tx, ty, objects, segmap):
    """
    Return the segmap label that belongs to the target at pixel position (tx, ty).
    0 means the target was not detected (background or outside map).
    """
    ny, nx = segmap.shape
    ix, iy = int(round(tx)), int(round(ty))

    # 1) direct pixel lookup
    if 0 <= iy < ny and 0 <= ix < nx:
        label = int(segmap[iy, ix])
        if label > 0:
            return label

    # 2) find nearest extracted source centroid within 5 pixels
    if len(objects) > 0:
        dists = np.hypot(objects["x"] - tx, objects["y"] - ty)
        j = int(np.argmin(dists))
        if dists[j] < 5.0:
            return j + 1   # sep labels are 1-indexed

    return 0   # not detected


def _contamination_from_segmap(bkg_sub, segmap, tx, ty, target_label,
                                aperture_radius_pix):
    """
    Pixel-level contamination fraction within a circular aperture.

    contamination = flux from non-target pixels inside aperture
                    ─────────────────────────────────────────────
                    total flux inside aperture  (target + others)

    Returns (contamination_fraction, total_aperture_flux).
    contamination_fraction is NaN when the aperture is empty or
    the total flux is non-positive.
    """
    ny, nx = bkg_sub.shape
    yy, xx = np.ogrid[:ny, :nx]
    in_aperture = (xx - tx) ** 2 + (yy - ty) ** 2 <= aperture_radius_pix ** 2

    total_flux = float(np.sum(bkg_sub[in_aperture]))

    if total_flux <= 0 or target_label == 0:
        return np.nan, total_flux

    # Only count flux from OTHER detected sources (label > 0), not
    # unlabeled background pixels whose noise would inflate contamination.
    contam_mask = in_aperture & (segmap > 0) & (segmap != target_label)
    contam_flux = float(np.sum(bkg_sub[contam_mask]))

    contamination = np.clip(contam_flux / total_flux, 0.0, 1.0)
    return contamination, total_flux


def _contamination_geometric(bkg_sub, objects, target_idx, tx, ty,
                              aperture_radius_pix, blend_sep_pix):
    """
    Fallback contamination estimate when no segmentation map is available.
    Uses the geometric overlap between each neighbor's circular aperture and
    the target's aperture to weight the neighbor's sep AUTO flux.

    contamination ≈ Σ_k (F_k × overlap_fraction(d_k, r))
                    ──────────────────────────────────────
                    F_target + Σ_k (F_k × overlap_fraction(d_k, r))
    """
    r = aperture_radius_pix
    if len(objects) == 0:
        return np.nan

    def _overlap_frac(d, r_):
        """Area fraction of circle(r_) at distance d that overlaps with circle(r_) at origin."""
        if d >= 2 * r_:
            return 0.0
        if d <= 0.0:
            return 1.0
        arg = np.clip(d / (2 * r_), -1, 1)
        A = 2 * r_ ** 2 * np.arccos(arg) - (d / 2) * np.sqrt(max(0, 4 * r_ ** 2 - d ** 2))
        return float(np.clip(A / (np.pi * r_ ** 2), 0.0, 1.0))

    f_target = float(objects["flux"][target_idx]) if target_idx >= 0 else 0.0

    contam_weighted = 0.0
    for k, obj in enumerate(objects):
        if k == target_idx:
            continue
        d = float(np.hypot(obj["x"] - tx, obj["y"] - ty))
        if d > blend_sep_pix:
            continue
        w = _overlap_frac(d, r)
        contam_weighted += float(obj["flux"]) * w

    denom = f_target + contam_weighted
    if denom <= 0:
        return np.nan
    return np.clip(contam_weighted / denom, 0.0, 1.0)


def _analyze_single_band(ra_i, dec_i, data_full, wcs_obj,
                          aperture_radius_arcsec, blend_sep_arcsec,
                          cutout_size_arcsec, detection_sigma, npixels,
                          deblend_nthresh, deblend_cont,
                          bkg_box_size, bkg_filtersize):
    """
    Run full blending analysis for one galaxy in one band.

    Returns
    -------
    contamination  : float   (NaN on failure)
    nearest_sep    : float   arcsec to nearest other source (np.inf if none)
    n_neighbors    : int     sources within blend_sep_arcsec
    is_blended     : bool
    """
    _nan = (np.nan, np.inf, 0, False)

    # ── pixel scale ──────────────────────────────────────────────────────────
    pix_scale          = _pixel_scale_arcsec(wcs_obj)          # arcsec/pix
    aperture_pix       = aperture_radius_arcsec / pix_scale
    blend_pix          = blend_sep_arcsec       / pix_scale
    cutout_pix         = int(np.ceil(cutout_size_arcsec / pix_scale))
    cutout_pix        += cutout_pix % 2          # make even

    # ── cutout ────────────────────────────────────────────────────────────────
    coord = SkyCoord(ra=ra_i * u.degree, dec=dec_i * u.degree)
    try:
        cutout = Cutout2D(data_full, coord,
                          size=cutout_pix,   # pixels
                          wcs=wcs_obj,
                          mode="partial",
                          fill_value=0.0)
    except (NoOverlapError, Exception) as exc:
        warnings.warn(f"Cutout failed for (RA={ra_i:.4f}, Dec={dec_i:.4f}): {exc}")
        return _nan

    stamp     = _sep_array(cutout.data)
    wcs_cut   = cutout.wcs

    # Target pixel coordinates in the cutout frame (x=col, y=row, scalar floats)
    pix = wcs_cut.world_to_pixel(coord)
    tx, ty = float(np.asarray(pix[0])), float(np.asarray(pix[1]))

    # ── background subtraction ────────────────────────────────────────────────
    box = bkg_box_size if bkg_box_size is not None else max(4, cutout_pix // 4)
    bkg_sub, rms = _background_subtract(stamp, box, bkg_filtersize)

    # ── source detection ──────────────────────────────────────────────────────
    objects, segmap = _extract_sources(bkg_sub, rms,
                                       detection_sigma, npixels,
                                       deblend_nthresh, deblend_cont)

    if len(objects) == 0:
        return _nan

    # ── match target to detected source ───────────────────────────────────────
    dists_pix = np.hypot(objects["x"] - tx, objects["y"] - ty)
    target_idx = int(np.argmin(dists_pix))
    # Accept if within blend_pix; otherwise target is undetected
    if dists_pix[target_idx] > blend_pix:
        target_idx = -1

    # ── nearest OTHER source ──────────────────────────────────────────────────
    other_dists_pix = dists_pix.copy()
    if target_idx >= 0:
        other_dists_pix[target_idx] = np.inf
    nearest_pix = float(np.min(other_dists_pix)) if len(objects) > 1 else np.inf
    nearest_sep = nearest_pix * pix_scale         # arcsec

    # ── neighbors within blend radius ─────────────────────────────────────────
    neighbor_mask = other_dists_pix <= blend_pix
    n_neighbors   = int(np.sum(neighbor_mask))
    is_blended    = n_neighbors > 0

    # ── contamination fraction ────────────────────────────────────────────────
    if segmap is not None:
        target_label  = _find_target_label(tx, ty, objects, segmap)
        contamination, _ = _contamination_from_segmap(
            bkg_sub, segmap, tx, ty, target_label, aperture_pix)
    else:
        # geometric fallback
        contamination = _contamination_geometric(
            bkg_sub, objects, target_idx, tx, ty, aperture_pix, blend_pix)

    return contamination, nearest_sep, n_neighbors, is_blended


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def check_blending(
    ra,
    dec,
    images,
    aperture_radius_arcsec=1.0,
    blend_sep_arcsec=2.0,
    cutout_size_arcsec=12.0,
    detection_sigma=2.0,
    npixels=5,
    deblend_nthresh=32,
    deblend_cont=0.005,
    multiband=True,
    bkg_box_size=None,
    bkg_filtersize=3,
):
    """
    Determine blending status and flux contamination for one or more galaxies.

    Parameters
    ----------
    ra : float or array-like
        Right ascension(s) of target galaxy(ies) in degrees.
    dec : float or array-like
        Declination(s) of target galaxy(ies) in degrees.
    images : tuple or dict
        Image data and WCS for one or more bands.

        Single band  →  (data_2d, wcs_object)
                        where data_2d is a 2-D numpy array and wcs_object is
                        an astropy.wcs.WCS instance.

        Multi-band   →  {"F200W": (data_2d, wcs_object),
                          "F444W": (data_2d, wcs_object), ...}

    aperture_radius_arcsec : float
        Radius (arcsec) of the circular aperture used for flux measurements.
        Should be comparable to the PSF FWHM (~0.06–0.16" for JWST NIRCam).
        Default: 1.0 arcsec.
    blend_sep_arcsec : float
        Neighbor search radius (arcsec).  Sources beyond this distance are
        ignored when counting blending neighbors.  Default: 2.0 arcsec.
    cutout_size_arcsec : float
        Side length (arcsec) of the postage stamp extracted around each target
        for source detection.  Should be large enough to capture all nearby
        neighbors.  Default: 12.0 arcsec.
    detection_sigma : float
        Source detection threshold in units of local background RMS.
        Lowering this finds fainter neighbors but increases false positives.
        Default: 2.0.
    npixels : int
        Minimum number of connected pixels required to define a source.
        Default: 5.
    deblend_nthresh : int
        Number of deblending sub-thresholds (sep / SExtractor parameter).
        Higher values produce finer deblending.  Default: 32.
    deblend_cont : float
        Minimum contrast ratio for deblending (sep / SExtractor parameter).
        Smaller values allow more aggressive deblending of close pairs.
        Default: 0.005.
    multiband : bool
        If True and `images` is a multi-band dict, run the analysis
        independently in each band and return per-band contamination.
        If False, only the first band is used.  Default: True.
    bkg_box_size : int or None
        Background mesh cell size in pixels.  If None, set automatically to
        ¼ of the cutout size.  Default: None.
    bkg_filtersize : int
        Size of the median filter applied to the background mesh (sep
        parameter).  Default: 3.

    Returns
    -------
    results : dict
        "is_blended" : np.ndarray of bool, shape (N,)
            True if the galaxy has ≥1 detected neighbor within
            `blend_sep_arcsec` in ANY analyzed band.

        "contamination_fraction" : np.ndarray of float
            Fraction of the flux inside the circular aperture that comes
            from pixels NOT belonging to the target galaxy.
            Defined as (total_aperture_flux – target_aperture_flux) /
                        total_aperture_flux,
            clipped to [0, 1].  NaN when the target is undetected.

            Shape: (N,) when a single band is analyzed.
                   (N, Nbands) when multiband=True and multiple bands
                   are provided.  Column order matches ``band_names``.

        "nearest_sep_arcsec" : np.ndarray of float, shape (N,)
            Angular separation (arcsec) to the nearest OTHER detected
            source.  np.inf if no other source is found.
            Minimum across bands when multiband=True.

        "n_blending_neighbors" : np.ndarray of int, shape (N,)
            Number of detected neighbors within `blend_sep_arcsec`.
            Maximum across bands when multiband=True.

        "band_names" : list of str
            Names of the analyzed bands, in the order they appear as
            columns of ``contamination_fraction`` (multiband case).

    Notes
    -----
    The pixel-level contamination calculation requires a segmentation map
    from sep.extract (available in sep ≥ 1.0).  If sep was compiled without
    this feature, a geometric circle-overlap approximation is used instead.

    All source detection is performed on the local, background-subtracted
    stamp.  The background is estimated with sep.Background using a mesh of
    size `bkg_box_size`.
    """
    # ── normalise coordinate inputs ───────────────────────────────────────────
    ra  = np.atleast_1d(np.asarray(ra,  dtype=float))
    dec = np.atleast_1d(np.asarray(dec, dtype=float))
    if ra.shape != dec.shape:
        raise ValueError("`ra` and `dec` must have the same shape.")
    N = len(ra)

    # ── normalise image inputs ────────────────────────────────────────────────
    if isinstance(images, tuple) and len(images) == 2:
        # Could be (data, wcs) OR ((data1, wcs1), (data2, wcs2))
        # Distinguish: first element should be a numpy-array-like for single band
        first = images[0]
        if isinstance(first, np.ndarray) or (hasattr(first, "shape") and hasattr(first, "ndim")):
            images = {"band": images}
        # else already a dict-like — handled below

    if not isinstance(images, dict):
        raise TypeError(
            "`images` must be a (data, wcs) tuple or a dict of {band: (data, wcs)}."
        )

    if not multiband or len(images) == 1:
        # use only the first band
        first_key = next(iter(images))
        images = {first_key: images[first_key]}

    band_names = list(images.keys())
    Nb = len(band_names)

    # ── output arrays ─────────────────────────────────────────────────────────
    contam_all   = np.full((N, Nb), np.nan)
    nearest_all  = np.full((N, Nb), np.inf)
    n_neigh_all  = np.zeros((N, Nb), dtype=int)

    # ── main loop ─────────────────────────────────────────────────────────────
    for b, band in enumerate(band_names):
        data_full, wcs_obj = images[band]

        for i in range(N):
            cont, near, n_nb, _ = _analyze_single_band(
                ra[i], dec[i], data_full, wcs_obj,
                aperture_radius_arcsec = aperture_radius_arcsec,
                blend_sep_arcsec       = blend_sep_arcsec,
                cutout_size_arcsec     = cutout_size_arcsec,
                detection_sigma        = detection_sigma,
                npixels                = npixels,
                deblend_nthresh        = deblend_nthresh,
                deblend_cont           = deblend_cont,
                bkg_box_size           = bkg_box_size,
                bkg_filtersize         = bkg_filtersize,
            )
            contam_all[i, b]  = cont
            nearest_all[i, b] = near
            n_neigh_all[i, b] = n_nb

    # ── aggregate across bands ────────────────────────────────────────────────
    is_blended        = np.any(n_neigh_all > 0, axis=1)          # (N,)
    nearest_sep       = np.min(nearest_all, axis=1)               # (N,)
    n_neighbors       = np.max(n_neigh_all, axis=1)               # (N,)

    if Nb == 1:
        contamination = contam_all[:, 0]                          # (N,)
    else:
        contamination = contam_all                                 # (N, Nb)

    return {
        "is_blended":             is_blended,
        "contamination_fraction": contamination,
        "nearest_sep_arcsec":     nearest_sep,
        "n_blending_neighbors":   n_neighbors,
        "band_names":             band_names,
    }
