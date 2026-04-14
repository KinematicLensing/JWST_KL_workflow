"""
background.py
-------------
NIRCam WFSS grism background processing:

1. ``background_grism_stage2``   – build and save median super-sky backgrounds.
2. ``assignwcs_grism_stage2``    – assign WCS + apply flat-field (→ lv1.5 files).
3. ``get_crds_dict_from_fits_header`` – helper to build a CRDS query dict.
4. ``my_grism_bkg_subtraction``  – per-exposure background subtraction including
   1/f noise removal and 2-D SExtractor background modelling.

Design note on global state
~~~~~~~~~~~~~~~~~~~~~~~~~~~
The original notebook stored ``cali_support_dir`` as a global variable.  Here
every function that needs it accepts it as an explicit parameter (or reads it
from a ``PipelineConfig`` instance).  The two processing functions
(``background_grism_stage2`` and ``my_grism_bkg_subtraction``) accept
``cali_support_dir`` directly so they can be used with ``multiprocessing.Pool``.
"""

from __future__ import annotations

import os
import time

import numpy as np
import matplotlib
matplotlib.use("Agg")           # non-interactive backend safe for pool workers
import matplotlib.pyplot as plt

from astropy.io import fits
from astropy.stats import sigma_clipped_stats, SigmaClip
from photutils.background import Background2D, SExtractorBackground
from photutils.segmentation import detect_sources

from nircam_wfss.noise import my_1overf_subtraction
from nircam_wfss.plotting import corner_text


# ---------------------------------------------------------------------------
# Stage-2 helpers
# ---------------------------------------------------------------------------

def get_crds_dict_from_fits_header(tmp_hd: fits.Header) -> dict:
    """
    Build a minimal CRDS query dictionary from a FITS primary header.

    Parameters
    ----------
    tmp_hd:
        FITS primary header of a NIRCam exposure.

    Returns
    -------
    crds_dict:
        Dictionary suitable for passing to ``crds.getreferences()``.
    """
    crds_dict: dict = {}
    crds_dict["INSTRUME"] = tmp_hd["INSTRUME"].upper()
    crds_dict["READPATT"] = tmp_hd["READPATT"].upper()
    crds_dict["SUBARRAY"] = tmp_hd["SUBARRAY"].upper()
    crds_dict["DATE-OBS"] = tmp_hd["DATE-OBS"]
    crds_dict["TIME-OBS"] = tmp_hd["TIME-OBS"]
    crds_dict["DETECTOR"] = tmp_hd["DETECTOR"]
    if crds_dict["INSTRUME"] == "NIRCAM":
        crds_dict["CHANNEL"] = (
            "LONG" if crds_dict["DETECTOR"] in ("NRCALONG", "NRCBLONG") else "SHORT"
        )
    crds_dict["MODULE"] = tmp_hd["MODULE"]
    crds_dict["EXP_TYPE"] = tmp_hd["EXP_TYPE"]
    crds_dict["FILTER"] = tmp_hd["FILTER"]
    crds_dict["PUPIL"] = tmp_hd["PUPIL"]
    return crds_dict


def background_grism_stage2(
    fanmes: np.ndarray,
    tmp_module: str,
    tmp_pupil: str,
    tmp_filter: str,
    cali_support_dir: str,
    overwrite: bool = False,
) -> None:
    """
    Create and save the sigma-clipped median super-sky background for one
    module / pupil / filter combination.

    After computing the raw sigma-clipped median, the background is post-processed
    to produce a smoother model:

    * GRISMR: the 2-D background is replaced by its 1-D column median, removing
      cross-dispersion row-to-row variation while preserving the detector-column
      illumination profile.
    * GRISMC: a transition-region (coronagraph boundary) is detected and handled
      separately for wavelengths below / above 4 µm, and a final (1, 31) median
      filter is applied along the dispersion axis.

    Parameters
    ----------
    arr_grism_rate:
        3-D stack of grism images, shape (ny, nx, n_frames).
    tmp_module:
        Detector module: ``'A'`` or ``'B'``.
    tmp_pupil:
        Grism pupil: ``'GRISMR'`` or ``'GRISMC'``.
    tmp_filter:
        NIRCam filter name (e.g. ``'F444W'``).
    cali_support_dir:
        Output directory for the background FITS file.
    overwrite:
        If ``False`` and the file already exists, skip recomputing.
    """
    from scipy import ndimage

    print(
        "\n>> Create Median Background for %s - module %s - grism %s"
        % (tmp_filter, tmp_module, tmp_pupil)
    )
    tmp_med_bkg_path = os.path.join(
        cali_support_dir,
        "median_bkg_%s_mod%s_%s.fits" % (tmp_filter, tmp_module, tmp_pupil),
    )
    if os.path.exists(tmp_med_bkg_path) and not overwrite:
        print(
            "Background model (%s-%s-%s) exists: %s"
            % (tmp_filter, tmp_module, tmp_pupil, tmp_med_bkg_path)
        )
        return
    
    arr_grism_rate = np.stack([
        fits.getdata(fanme) for fanme in fanmes
    ], axis=-1)

    # 2-sigma-clipped median along the frame axis
    tmp_med_bkg = sigma_clipped_stats(
        arr_grism_rate, axis=-1, sigma_upper=2.0, maxiters=10
    )[1]

    # Post-process the median background to produce a smooth model
    tmp_grism_bkg = tmp_med_bkg.copy()

    if tmp_pupil[-1].upper() == "R":
        # GRISMR disperses horizontally: replace with 1-D column median to
        # remove cross-dispersion row variation, keeping the column profile.
        tmp_1d_med_bkg = np.nanmedian(tmp_grism_bkg, axis=0)
        tmp_grism_bkg = np.ones_like(tmp_grism_bkg) * tmp_1d_med_bkg

    elif tmp_pupil[-1].upper() == "C":
        # GRISMC disperses vertically: handle the coronagraph transition region.
        tmp_grism_bkg_new = tmp_grism_bkg.copy()

        # Detect the transition row (where the background level stops being flat)
        y_med, y_std = sigma_clipped_stats(tmp_grism_bkg, axis=1)[1:]
        y_trans = int(np.max(np.where(y_med / np.nanpercentile(y_med, 95) > 0.99))) - 100

        # Hard-coded overrides for known filter-specific transition positions
        if tmp_filter == "F322W2":
            y_trans = 750
        elif tmp_filter == "F444W":
            y_trans = 1000

        print("   GRISMC transition row: y_trans = %d" % y_trans)

        if int(tmp_filter[1:4]) < 400:
            # Short-wavelength filters (e.g. F322W2): smooth region is below y_trans
            tmp_grism_bkg_new[:y_trans] = (np.ones_like(tmp_grism_bkg) * y_med).T[:y_trans]
            # Coronagraph transmission region above y_trans
            tmp_grism_bkg_new_top = ndimage.median_filter(
                np.nan_to_num(tmp_grism_bkg_new[y_trans:]), (11, 1), mode="constant"
            ).copy()
            scale_top = (
                sigma_clipped_stats(tmp_grism_bkg_new_top[0:15], axis=0)[1]
                / np.nanmedian(tmp_grism_bkg_new[y_trans - 10 : y_trans - 5], axis=0)
            )
            tmp_grism_bkg_new[y_trans:] = tmp_grism_bkg_new_top / scale_top
        else:
            # Long-wavelength filters (e.g. F444W): smooth region is above y_trans
            tmp_grism_bkg_new[y_trans:] = (np.ones_like(tmp_grism_bkg) * y_med).T[y_trans:]
            # Coronagraph transmission region below y_trans
            tmp_grism_bkg_new_bottom = ndimage.median_filter(
                np.nan_to_num(tmp_grism_bkg_new[:y_trans]), (11, 1), mode="constant"
            ).copy()
            scale_top = (
                sigma_clipped_stats(tmp_grism_bkg_new_bottom[-15:], axis=0)[1]
                / np.nanmedian(tmp_grism_bkg_new[y_trans + 5 : y_trans + 10], axis=0)
            )
            tmp_grism_bkg_new[:y_trans] = tmp_grism_bkg_new_bottom / scale_top

        # Final smooth along the dispersion axis
        tmp_grism_bkg = ndimage.median_filter(tmp_grism_bkg_new, size=(1, 31))

    hdul = fits.HDUList(fits.PrimaryHDU(tmp_grism_bkg))
    hdul[0].header["OBJECT"] = "BKG_%s_mod%s_%s" % (tmp_filter, tmp_module, tmp_pupil)
    hdul[0].header["FILTER"] = tmp_filter
    hdul[0].header["MODULE"] = tmp_module
    hdul[0].header["PUPIL"]  = tmp_pupil
    hdul.writeto(tmp_med_bkg_path, overwrite=True)
    print(">> Saved background model to %s" % tmp_med_bkg_path)


def assignwcs_grism_stage2(
    rate_grism_file: str,
    v1p5_dir: str,
    overwrite: bool = False,
) -> int:
    """
    Assign WCS and apply flat-field to a NIRCam grism *rate* file, producing
    a level-1.5 FITS file.

    The function "tricks" the JWST pipeline into treating the grism exposure
    as a direct image so that the standard ``AssignWcsStep`` and flat-field
    step can be applied without grism-throughput contamination.

    Parameters
    ----------
    rate_grism_file:
        Path to a NIRCam grism *rate.fits* file (stage-1 output).
    v1p5_dir:
        Output directory for v1.5 calibrated files.
        The output file is named ``<basename>_rate_lv1.5.fits``.
    overwrite:
        If ``False`` and the output file exists, skip processing.

    Returns
    -------
    status:
        0 if skipped (file exists), 1 if processed successfully.
    """
    import crds
    from jwst import assign_wcs, datamodels
    from jwst.flatfield import flat_field

    grism_save_path = os.path.join(
        v1p5_dir,
        os.path.basename(rate_grism_file).replace("rate.fits", "rate_lv1.5.fits"),
    )
    if os.path.exists(grism_save_path) and not overwrite:
        print("WCS-assigned file exists: %s" % grism_save_path)
        return 0

    with fits.open(rate_grism_file) as rate_grism_fits:
        tmp_grism_hd = rate_grism_fits[0].header
        tmp_filter = tmp_grism_hd["FILTER"]
        tmp_module = tmp_grism_hd["MODULE"]
        tmp_pupil  = tmp_grism_hd["PUPIL"]

        print(
            "   read %s (%s - %s - %s)"
            % (rate_grism_file, tmp_filter, tmp_module, tmp_pupil)
        )

        # Assign WCS (fool pipeline by labelling as NRC_IMAGE)
        siaf_file = crds.getreferences(
            tmp_grism_hd, reftypes=["distortion"], ignore_cache=False
        )["distortion"]
        grism_wcs_step = assign_wcs.assign_wcs_step.AssignWcsStep(
            override_distortion=siaf_file
        )
        rate_grism_fits[0].header["EXP_TYPE"] = "NRC_IMAGE"
        grism_image = datamodels.image.ImageModel(rate_grism_fits)
        grism_with_wcs = grism_wcs_step.run(grism_image)

        # Flat-field (use CLEAR pupil to avoid grism throughput)
        date_obs = tmp_grism_hd["DATE"].split("T")[0]
        time_obs = tmp_grism_hd["DATE"].split("T")[1]
        tmp_flat_path = crds.getreferences(
            {
                "INSTRUME": tmp_grism_hd["INSTRUME"],
                "READPATT": tmp_grism_hd["READPATT"],
                "SUBARRAY": tmp_grism_hd["SUBARRAY"],
                "DATE-OBS": date_obs,
                "TIME-OBS": time_obs,
                "DETECTOR": tmp_grism_hd["DETECTOR"],
                "CHANNEL":  tmp_grism_hd["CHANNEL"],
                "MODULE":   tmp_grism_hd["MODULE"],
                "EXP_TYPE": "NRC_IMAGE",
                "FILTER":   tmp_grism_hd["FILTER"],
                "PUPIL":    "CLEAR",
            }
        )["flat"]
        try:
            flat_field.do_flat_field(grism_with_wcs, datamodels.FlatModel(tmp_flat_path))
        except Exception as e:
            print("Error applying flat-field %s to %s: %s" % (tmp_flat_path, rate_grism_file, str(e)))
        print("Saving lv1.5 file to:", grism_save_path)
        grism_with_wcs.save(grism_save_path)

    return 1


# ---------------------------------------------------------------------------
# Per-exposure background subtraction
# ---------------------------------------------------------------------------

def my_grism_bkg_subtraction(
    rate_grism_file: str,
    cali_support_dir: str,
    plot_dir: str | None = None,
    use_robust_background: bool = True,
) -> None:
    """
    Apply background subtraction to one lv1.5 grism exposure file, in place.

    Processing steps:
    1. Load the median super-sky for the matching filter/module/pupil.
    2. Scale and subtract the super-sky.
    3. Run source detection and estimate the residual 2-D SExtractor background.
    4. Subtract the residual 2-D background.
    5. Remove 1/f noise stripes (per-row / per-column median).
    6. Overwrite the SCI extension in the original file.

    Parameters
    ----------
    rate_grism_file:
        Path to the lv1.5 grism FITS file.  Modified in place.
    cali_support_dir:
        Directory where calibration support median super-sky FITS files are stored.
    plot_dir:
        If provided, diagnostic PNG plots are saved here.  If ``None``,
        no plots are produced.
    use_robust_background:
        If ``True``, apply an additional robust median background subtraction 
        after 1/f noise removal. Otherwise, subtract background estimated from 
        2-D SExtractor before 1/f noise removal. 
        This can help mitigate any residual background offsets.
    """
    print(">>> subtract background for %s" % os.path.basename(rate_grism_file))

    tmp_grism_fits = fits.open(rate_grism_file, mode="update")
    tmp_grism_img  = tmp_grism_fits[1].data          # SCI frame
    tmp_grism_hd   = tmp_grism_fits[0].header

    tmp_filter = tmp_grism_hd["FILTER"]
    tmp_module = tmp_grism_hd["MODULE"]
    tmp_pupil  = tmp_grism_hd["PUPIL"]

    # --- Step 1: load and scale median super-sky ---
    bkg_path = os.path.join(
        cali_support_dir,
        "median_bkg_%s_mod%s_%s.fits" % (tmp_filter, tmp_module, tmp_pupil),
    )
    tmp_grism_bkg = fits.getdata(bkg_path)

    if tmp_pupil == "GRISMR":
        tmp_1d_med_bkg = np.nanmedian(tmp_grism_bkg, axis=0)
        tmp_1d_sci_bkg = sigma_clipped_stats(
            tmp_grism_img, axis=0, sigma_upper=2.5
        )[1]
        tmp_scale_bkg = np.nanmedian(
            (tmp_1d_sci_bkg / tmp_1d_med_bkg)[np.isfinite(tmp_1d_med_bkg)]
        )
        tmp_grism_bkg = tmp_grism_bkg * tmp_scale_bkg

    elif tmp_pupil == "GRISMC":
        arg_scale = (
            np.arange(600, 2048)
            if tmp_filter == "F480M"
            else np.arange(2048)
        )
        tmp_1d_med_bkg = np.nanmedian(tmp_grism_bkg[arg_scale], axis=1)
        tmp_1d_sci_bkg = sigma_clipped_stats(
            tmp_grism_img[arg_scale], axis=1, sigma_upper=2.5
        )[1]
        tmp_scale_bkg = np.nanmedian(
            (tmp_1d_sci_bkg / tmp_1d_med_bkg)[np.isfinite(tmp_1d_med_bkg)]
        )
        tmp_grism_bkg = tmp_grism_bkg * tmp_scale_bkg

    # --- Optional diagnostic plot (before/after state) ---
    if plot_dir is not None:
        plt.close("all")
        fig, ax = plt.subplots(2, 2, figsize=(15, 15))
        ax = ax.flatten()
        ax[0].imshow(
            np.nan_to_num(tmp_grism_img), origin="lower",
            vmin=np.nanpercentile(tmp_grism_img, 2.5),
            vmax=np.nanpercentile(tmp_grism_img, 97.5),
        )
        corner_text(ax[0], s="(1) Flat-Fielded", color="w", loc=2,
                    weight="bold", fontsize=20)

    # --- Step 2: subtract scaled super-sky ---
    tmp_grism_img = tmp_grism_img - tmp_grism_bkg
    tmp_grism_fits[0].header["HISTORY"] = (
        "sigma-clipped median sky background subtracted on %s"
        % time.strftime("%Y/%m/%d", time.localtime())
    )
    if plot_dir is not None:
        ax[1].imshow(
            np.nan_to_num(tmp_grism_img), origin="lower",
            vmin=np.nanpercentile(tmp_grism_img, 2.5),
            vmax=np.nanpercentile(tmp_grism_img, 97.5),
        )
        corner_text(ax[1], s="(2) After super-sky subtraction", color="w", loc=2,
                    weight="bold", fontsize=20)

    if not use_robust_background:
        # --- Step 3: detect sources; residual 2-D background ---
        _, tmp_med, tmp_rms = sigma_clipped_stats(tmp_grism_img[100:700, 100:700])
        segment_map = detect_sources(tmp_grism_img - tmp_med, tmp_rms * 1.0, npixels=100)
        sigma_clip    = SigmaClip(sigma=2.5)
        bkg_estimator = SExtractorBackground()
        #tmp_mask = np.zeros_like(segment_map.data, dtype=bool)
        bkg = Background2D(
            tmp_grism_img, (64, 64),
            filter_size=(9, 9),
            mask=segment_map,
            sigma_clip=sigma_clip,
            bkg_estimator=bkg_estimator,
        )
        # --- Step 4: subtract residual 2-D background ---
        tmp_grism_img = tmp_grism_img - bkg.background
        if plot_dir is not None:
            ax[2].imshow(
                np.nan_to_num(tmp_grism_img), origin="lower",
                vmin=np.nanpercentile(tmp_grism_img, 2.5),
                vmax=np.nanpercentile(tmp_grism_img, 97.5),
            )
            corner_text(ax[2], s="(3) After 2D background subtraction",
                        color="w", loc=2, weight="bold", fontsize=20)

    # --- Step 5: 1/f noise removal ---
    _, tmp_med, tmp_rms = sigma_clipped_stats(tmp_grism_img[100:1700, 100:1700].flatten()[::10])
    segment_map_2 = detect_sources(tmp_grism_img, tmp_rms * 2.0, npixels=200)
    backg_mask_1f = np.ones_like(tmp_grism_img, dtype=int)
    if segment_map_2 is not None:
        backg_mask_1f[segment_map_2.data > 0] = 0

    if tmp_pupil[-1] == "C":
        # Subtract column-wise stripes (GRISMC disperses vertically)
        tmp_grism_img, _ = my_1overf_subtraction(
            tmp_grism_img, backg_mask=backg_mask_1f, amplifiers=4
        )
    elif tmp_pupil[-1] == "R":
        # Subtract row-wise stripes (GRISMR disperses horizontally)
        tmp_col_sub, _ = my_1overf_subtraction(
            tmp_grism_img.T, backg_mask=backg_mask_1f.T, amplifiers=1
        )
        tmp_grism_img = tmp_col_sub.T
    if plot_dir is not None:
        ax[3].imshow(
            np.nan_to_num(tmp_grism_img), origin="lower",
            vmin=np.nanpercentile(tmp_grism_img, 2.5),
            vmax=np.nanpercentile(tmp_grism_img, 97.5),
        )
        corner_text(ax[3], s="(4) After 1/f subtraction",
                    color="w", loc=2, weight="bold", fontsize=20)

    if use_robust_background:
        tmp_grism_img = tmp_grism_img - robust_median_bkg(tmp_grism_img)
        if plot_dir is not None:
            ax[2].imshow(
                np.nan_to_num(tmp_grism_img), origin="lower",
                vmin=np.nanpercentile(tmp_grism_img, 2.5),
                vmax=np.nanpercentile(tmp_grism_img, 97.5),
            )
            corner_text(ax[2], s="(3) After 2D background subtraction",
                        color="w", loc=2, weight="bold", fontsize=20)

    # --- Step 6: write back to file ---
    tmp_grism_fits[1].data = tmp_grism_img.astype(np.float32)
    tmp_grism_fits[0].header["HISTORY"] = (
        "1/f noise + SExtractor 2D background subtracted on %s"
        % time.strftime("%Y/%m/%d", time.localtime())
    )

    # --- Optional: save diagnostic plot ---
    if plot_dir is not None:
        plt.tight_layout()
        out_plot = os.path.join(
            plot_dir,
            "bkg_%s" % os.path.basename(rate_grism_file).replace(".fits", ".png"),
        )
        fig.savefig(out_plot, dpi=80, bbox_inches="tight")
        plt.close("all")

    tmp_grism_fits.flush()
    tmp_grism_fits.close()


def robust_median_bkg(img, mask = False):
    '''get robust background of image by masking brightest 25th percentile of the data in row / column '''
    arr_x_avg = np.nanpercentile(img, 75, axis = 0)
    arr_y_avg = np.nanpercentile(img, 75, axis = 1)
    img_mask = np.ones_like(img)
    img_mask[np.where(arr_y_avg < np.nanmedian(arr_y_avg))[0],:]  = 0
    img_mask[:,np.where(arr_x_avg < np.nanmedian(arr_x_avg))[0]] = 0
    tmp_global_bkg = sigma_clipped_stats(img[img_mask==0].flatten(), sigma_upper = 2., maxiters = 5)[1]
    if mask == True: return tmp_global_bkg, img_mask
    else: return tmp_global_bkg
