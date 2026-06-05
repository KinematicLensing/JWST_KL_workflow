"""
wcs_utils.py
------------
WCS and attitude-matrix utilities for JWST NIRCam astrometric corrections.

These functions translate between V2/V3 focal-plane coordinates and sky
coordinates (RA/Dec/PA), following JWST-STScI-001550 SM-12 §6.1.
"""

from __future__ import annotations

import numpy as np


# ---------------------------------------------------------------------------
# Fundamental rotation matrices
# ---------------------------------------------------------------------------

def rotate(axis: int, angle: float) -> np.ndarray:
    """
    Build a 3×3 rotation matrix about the given axis.

    Parameters
    ----------
    axis:
        Rotation axis: 1 (X), 2 (Y), or 3 (Z).
    angle:
        Rotation angle in degrees.

    Returns
    -------
    r:
        3×3 rotation matrix.
    """
    if axis not in (1, 2, 3):
        raise ValueError("axis must be 1, 2, or 3")
    r = np.zeros((3, 3))
    ax0 = axis - 1
    theta = np.deg2rad(angle)
    r[ax0, ax0] = 1.0
    ax1 = (ax0 + 1) % 3
    ax2 = (ax0 + 2) % 3
    r[ax1, ax1] =  np.cos(theta)
    r[ax2, ax2] =  np.cos(theta)
    r[ax1, ax2] = -np.sin(theta)
    r[ax2, ax1] =  np.sin(theta)
    return r


def attitude(
    v2: float,
    v3: float,
    ra: float,
    dec: float,
    pa: float,
) -> np.ndarray:
    """
    Build the attitude rotation matrix that maps V2/V3 to RA/Dec/PA.

    Rotates a unit vector at (v2, v3) to the pointing (ra, dec, pa)
    following JWST-STScI-001550 SM-12 §6.1.

    Parameters
    ----------
    v2, v3:
        Focal-plane V2 and V3 coordinates in arcseconds.
    ra, dec, pa:
        Pointing RA, Dec and position angle in degrees.

    Returns
    -------
    m:
        3×3 combined rotation matrix  M = Mra · Mdec · Mpa · Mv3 · Mv2.
    """
    v2d = v2 / 3600.0
    v3d = v3 / 3600.0

    mv2 = rotate(3, -v2d)
    mv3 = rotate(2,  v3d)
    mra = rotate(3,  ra)
    mdec = rotate(2, -dec)
    mpa = rotate(1,  -pa)

    m = np.dot(mv3, mv2)
    m = np.dot(mpa, m)
    m = np.dot(mdec, m)
    m = np.dot(mra, m)
    return m


# ---------------------------------------------------------------------------
# Astrometric offset correction helpers
# ---------------------------------------------------------------------------

def apply_astrometric_correction(
    grism_hd_sci: "astropy.io.fits.Header",  # noqa: F821
    sw_astrom_row: "astropy.table.Row",       # noqa: F821
) -> "astropy.io.fits.Header":               # noqa: F821
    """
    Apply RA/Dec shift and rotation correction derived from SW astrometry to
    a grism SCI header (modifying it in place).

    The correction updates CRVAL1, CRVAL2, and the CD matrix to account for
    the residual pointing offset measured from simultaneous SW direct images.

    Parameters
    ----------
    grism_hd_sci:
        Fits SCI-extension header of the grism exposure.  Modified in place.
    sw_astrom_row:
        A row from the SW astrometry table with columns
        ``dRA`` (arcsec), ``dDEC`` (arcsec), and ``theta`` (degrees).

    Returns
    -------
    grism_hd_sci:
        The same header object, modified in place.
    """
    crval2 = grism_hd_sci["CRVAL2"]

    # Shift CRVAL
    grism_hd_sci["CRVAL1"] -= (
        sw_astrom_row["dRA"] / np.cos(np.deg2rad(crval2)) / 3600.0
    )
    grism_hd_sci["CRVAL2"] -= sw_astrom_row["dDEC"] / 3600.0

    # Rotate CD matrix
    cd = np.array(
        [
            [grism_hd_sci["CD1_1"], grism_hd_sci["CD1_2"]],
            [grism_hd_sci["CD2_1"], grism_hd_sci["CD2_2"]],
        ]
    )
    theta = -sw_astrom_row["theta"]
    rot = np.array(
        [
            [np.cos(np.deg2rad(theta)), -np.sin(np.deg2rad(theta))],
            [np.sin(np.deg2rad(theta)),  np.cos(np.deg2rad(theta))],
        ]
    )
    cd = np.matmul(cd, rot)
    grism_hd_sci["CD1_1"] = cd[0, 0]
    grism_hd_sci["CD1_2"] = cd[0, 1]
    grism_hd_sci["CD2_1"] = cd[1, 0]
    grism_hd_sci["CD2_2"] = cd[1, 1]

    return grism_hd_sci
