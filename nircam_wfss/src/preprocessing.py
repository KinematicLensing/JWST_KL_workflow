"""
preprocessing.py
----------------
Per-exposure preprocessing steps applied to lv1.5 grism FITS files:

1. ``grism_hot_pix_rejection`` – flag hot pixels via a roll-difference test.
2. ``my_grism_cont_subtraction`` – separate continuum and emission-line images
   using a per-row / per-column median filter with a central hole.

Both functions modify the input FITS file in place and optionally save
diagnostic PNG plots.
"""

from __future__ import annotations

import os
import time

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from astropy.io import fits
from astropy.stats import sigma_clipped_stats

from scipy import ndimage

from nircam_wfss.noise import my_1overf_subtraction
from nircam_wfss.plotting import corner_text


# ---------------------------------------------------------------------------
# Hot pixel rejection
# ---------------------------------------------------------------------------

def grism_hot_pix_rejection(
    rate_grism_file: str,
    sigma_hot: float = 20.0,
    plot_dir: str | None = None,
) -> tuple[int, int]:
    """
    Identify and flag hot pixels in a lv1.5 grism FITS file.

    Uses a roll-difference method: a pixel is hot if the difference between
    it and its neighbour along the dispersion axis exceeds ``sigma_hot``
    times the robust standard deviation.  The SCI, ERR, and DQ extensions
    are updated in place.

    Parameters
    ----------
    rate_grism_file:
        Path to the lv1.5 grism FITS file.  Modified in place.
    sigma_hot:
        Detection threshold in units of the roll-difference standard
        deviation.  Default 20 (conservative).
    plot_dir:
        If provided, a diagnostic before/after PNG is saved here.

    Returns
    -------
    n_hot1, n_hot2:
        Number of hot pixels flagged in the first and second pass.
    """
    tmp_grism_fits = fits.open(rate_grism_file, mode="update")
    tmp_grism_img  = tmp_grism_fits[1].data.copy()
    tmp_grism_hd   = tmp_grism_fits[0].header

    # Roll direction depends on dispersion orientation
    tmp_axis = 1 if tmp_grism_hd["PUPIL"].upper() == "GRISMR" else 0

    # --- First pass: forward roll ---
    diff1 = np.nan_to_num(tmp_grism_img) - np.roll(
        np.nan_to_num(tmp_grism_img), 1, axis=tmp_axis
    )
    std_roll = sigma_clipped_stats(diff1.flatten(), sigma=5)[-1]
    y_hot_1, x_hot_1 = np.where(diff1 > sigma_hot * std_roll)

    tmp_grism_img_new = tmp_grism_img.copy()
    tmp_grism_img_new[y_hot_1, x_hot_1] = np.nan
    tmp_grism_fits["SCI"].data[y_hot_1, x_hot_1] = np.nan
    tmp_grism_fits["ERR"].data[y_hot_1, x_hot_1] = np.nan
    tmp_grism_fits["DQ"].data[y_hot_1, x_hot_1]  = 1

    # --- Second pass: backward roll ---
    diff2 = np.nan_to_num(tmp_grism_img_new) - np.roll(
        np.nan_to_num(tmp_grism_img_new), -1, axis=tmp_axis
    )
    y_hot_2, x_hot_2 = np.where(diff2 > sigma_hot * std_roll)

    tmp_grism_img_new[y_hot_2, x_hot_2] = np.nan
    tmp_grism_fits["SCI"].data[y_hot_2, x_hot_2] = np.nan
    tmp_grism_fits["ERR"].data[y_hot_2, x_hot_2] = np.nan
    tmp_grism_fits["DQ"].data[y_hot_2, x_hot_2]  = 1

    n_hot1, n_hot2 = len(y_hot_1), len(y_hot_2)
    tmp_grism_fits[0].header["HISTORY"] = (
        "%d additional hot pixels masked on %s"
        % (n_hot1 + n_hot2, time.strftime("%Y/%m/%d", time.localtime()))
    )
    print(
        "%s: %d + %d hot pixels flagged in 1st / 2nd pass"
        % (os.path.basename(rate_grism_file), n_hot1, n_hot2)
    )

    # --- Optional diagnostic plot ---
    if plot_dir is not None and (n_hot1 + n_hot2) > 0:
        idx_hot = np.random.choice(n_hot1 + n_hot2)
        all_x = np.concatenate((x_hot_1, x_hot_2))
        all_y = np.concatenate((y_hot_1, y_hot_2))
        tmp_cx, tmp_cy = all_x[idx_hot], all_y[idx_hot]

        fig, axes = plt.subplots(1, 2, figsize=(10, 5.5))
        for ax_i, img_i, title_i in zip(
            axes,
            [tmp_grism_img, tmp_grism_img_new],
            ["Before Masking", "After Masking"],
        ):
            ax_i.imshow(
                np.nan_to_num(img_i), origin="lower",
                vmin=-std_roll, vmax=std_roll * 3,
            )
            corner_text(ax_i, loc=1, s=title_i, color="w", weight="bold")
            ax_i.set_xticks([])
            ax_i.set_yticks([])

        plt.tight_layout()
        out_path = os.path.join(
            plot_dir,
            "hotpix_%s" % os.path.basename(rate_grism_file).replace(".fits", ".png"),
        )
        fig.savefig(out_path, dpi=80, bbox_inches="tight")
        plt.close("all")

    tmp_grism_fits.flush()
    tmp_grism_fits.close()
    return n_hot1, n_hot2


# ---------------------------------------------------------------------------
# Continuum subtraction / emission-line image creation
# ---------------------------------------------------------------------------

def my_grism_cont_subtraction(
    rate_grism_file: str,
    plot_dir: str | None = None,
) -> int:
    """
    Subtract the continuum from a background-subtracted lv1.5 grism exposure
    and append emission-line (``EMLINE``) and continuum (``CONT``) extensions.

    The algorithm uses a 2-pass approach:

    * Pass 1: median filter with a central hole along the dispersion direction
      to isolate emission-line flux, then remove 1/f residuals.
    * Pass 2: re-estimate the continuum after masking high-S/N emission-line
      pixels from Pass 1, producing a cleaner emission-line image.

    If the file already has the EMLINE extension (9 HDUs), the function
    returns immediately without reprocessing.

    Parameters
    ----------
    rate_grism_file:
        Path to the lv1.5 grism FITS file.  If the EMLINE extension is
        absent, it is appended in place.
    plot_dir:
        If provided, a 4-panel diagnostic PNG is saved here.

    Returns
    -------
    status:
        0 if the file already had the EMLINE extension (skipped),
        1 if processed successfully.
    """
    tmp_grism_fits = fits.open(rate_grism_file, mode="update")
    if len(tmp_grism_fits) >= 9:            # already processed (EMLINE ± CONT appended)
        tmp_grism_fits.close()
        return 0

    tmp_grism_sci = tmp_grism_fits["SCI"]
    tmp_grism_img = tmp_grism_sci.data[4:-4, 4:-4].astype(float)
    tmp_grism_hd  = tmp_grism_fits[0].header
    tmp_pupil_last = tmp_grism_hd["PUPIL"][-1]   # 'R' or 'C'

    # --- Build median-filter footprint (horizontal hole for R, vertical for C) ---
    L_box, L_mask = 25, 4
    if tmp_pupil_last == "R":
        mf_footprint = np.ones((1, L_box * 2 + 1))
        mf_footprint[:, L_box - L_mask : L_box + L_mask + 1] = 0
        pass2a_fp = np.ones((1, 150))
        pass2b_fp = np.ones((1, 51))
    elif tmp_pupil_last == "C":
        mf_footprint = np.ones((L_box * 2 + 1, 1))
        mf_footprint[L_box - L_mask : L_box + L_mask + 1, :] = 0
        pass2a_fp = np.ones((150, 1))
        pass2b_fp = np.ones((51, 1))
    else:
        raise KeyError("PUPIL must end in 'R' or 'C', got: %s" % tmp_grism_hd["PUPIL"])

    # --- Pass 1: estimate continuum via median filter ---
    tmp_cont = ndimage.median_filter(
        np.nan_to_num(tmp_grism_img), footprint=mf_footprint, mode="reflect"
    )
    tmp_emline_1 = tmp_grism_img - tmp_cont

    # 1/f removal on pass-1 emission-line image
    rowSub_h, _ = my_1overf_subtraction(tmp_emline_1, amplifiers=1)
    rowSub_v, _ = my_1overf_subtraction(rowSub_h.T, amplifiers=1)
    tmp_emline_pass1 = rowSub_v.T

    # --- Pass 2: mask bright emission pixels and re-estimate continuum ---
    tmp_rms = sigma_clipped_stats(tmp_emline_pass1)[-1]
    high_sn_mask = np.abs(tmp_emline_pass1 / tmp_rms) > 2.0

    arg_x_valid = np.where(
        np.sum(np.isnan(tmp_grism_img), axis=0) != tmp_grism_img.shape[0]
    )[0]

    tmp_img_copy = tmp_grism_img.copy()
    tmp_img_medflt_a = tmp_img_copy.copy()
    tmp_img_medflt_a[:, arg_x_valid] = ndimage.median_filter(
        tmp_img_copy[:, arg_x_valid], footprint=pass2a_fp, mode="reflect"
    )
    # Replace high-S/N emission pixels with the coarse continuum estimate
    tmp_img_copy[high_sn_mask] = tmp_img_medflt_a[high_sn_mask]

    tmp_img_medflt_b = tmp_img_medflt_a.copy()
    tmp_img_medflt_b[:, arg_x_valid] = ndimage.median_filter(
        tmp_img_copy[:, arg_x_valid], footprint=pass2b_fp, mode="reflect"
    )

    tmp_emline_2 = tmp_grism_img - np.nan_to_num(tmp_img_medflt_b)
    rowSub2_h, _ = my_1overf_subtraction(tmp_emline_2, amplifiers=1)
    rowSub2_v, _ = my_1overf_subtraction(rowSub2_h.T, amplifiers=1)
    tmp_emline_pass2 = rowSub2_v.T
    tmp_cont_final   = tmp_grism_img - tmp_emline_pass2

    # --- Optional diagnostic plot ---
    if plot_dir is not None:
        fig, ax = plt.subplots(2, 2, figsize=(14.5, 14.5))
        ax = ax.flatten()
        panels = [
            (tmp_grism_img,    "Original (BKG Subtracted)"),
            (tmp_cont_final,   "2nd-Pass Continuum Only"),
            (tmp_emline_pass1, "1st-Pass EMLINE Only"),
            (tmp_emline_pass2, "2nd-Pass EMLINE Only"),
        ]
        for tmp_ax, (img_i, label_i) in zip(ax, panels):
            tmp_ax.imshow(
                np.nan_to_num(img_i), origin="lower",
                vmin=np.nanpercentile(img_i, 5),
                vmax=np.nanpercentile(img_i, 95),
            )
            corner_text(tmp_ax, label_i, color="w", loc=2,
                        weight="bold", fontsize=20)
            tmp_ax.set_xticks([])
            tmp_ax.set_yticks([])
        if len(panels) > 2:
            corner_text(ax[2], "RMS=%.4f" % sigma_clipped_stats(tmp_emline_pass1)[-1],
                        color="w", loc=4, weight="bold", fontsize=20)
            corner_text(ax[3], "RMS=%.4f" % sigma_clipped_stats(tmp_emline_pass2)[-1],
                        color="w", loc=4, weight="bold", fontsize=20)
        plt.tight_layout()
        out_path = os.path.join(
            plot_dir,
            "cont_%s" % os.path.basename(rate_grism_file).replace(".fits", ".png"),
        )
        fig.savefig(out_path, dpi=80, bbox_inches="tight")
        plt.close("all")

    # --- Append EMLINE and CONT extensions ---
    # Pad back to full-frame size
    ny_full, nx_full = tmp_grism_fits["SCI"].data.shape
    emline_full = np.zeros((ny_full, nx_full), dtype=np.float32)
    cont_full   = np.zeros((ny_full, nx_full), dtype=np.float32)
    emline_full[4:-4, 4:-4] = tmp_emline_pass2
    cont_full[4:-4, 4:-4]   = tmp_cont_final

    hdu_emline = fits.ImageHDU(emline_full, name="EMLINE")
    hdu_emline.header["HISTORY"] = (
        "Continuum-subtracted emission-line image, created on %s"
        % time.strftime("%Y/%m/%d", time.localtime())
    )
    hdu_cont = fits.ImageHDU(cont_full, name="CONT")
    hdu_cont.header["HISTORY"] = (
        "Continuum image, created on %s"
        % time.strftime("%Y/%m/%d", time.localtime())
    )

    tmp_grism_fits.append(hdu_emline)
    tmp_grism_fits.append(hdu_cont)
    tmp_grism_fits.flush()
    tmp_grism_fits.close()
    return 1
