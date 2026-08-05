#!/usr/bin/env python
"""
run_pipeline.py
===============
Stand-alone NIRCam WFSS grism extraction pipeline script.

All adjustable parameters live in a YAML configuration file passed on the
command line.  A reference YAML file is provided in the package's
``configs/`` directory.

Pipeline stages
---------------
1.  Stage-2a   – assign WCS + flat-field grism rate files  → lv1.5 files.
2.  Stage-2b   – build median super-sky backgrounds.
3.  Stage-2c   – background subtraction + 1/f noise removal (per exposure).
4.  Stage-2d   – hot-pixel rejection.
5.  Stage-2e   – continuum subtraction → EMLINE extensions.
6.  Stage-3    – reduce SW direct images through calwebb_image2.
7.  Stage-4    – astrometry calibration / load astrometry table.
8.  Stage-5    – load source catalogue, build per-frame POM-applied catalogs.
9.  Stage-6    – extract 2-D grism spectra for every source.
10. Stage-7    – extract 1-D spectra + diagnostic plots.
11. Stage-8    – extract 2-D emission-line cutouts (requires z_grism).

Usage
-----
    python run_pipeline.py <path/to/config.yaml>

A reference configuration file is at::

    configs/PID3577_CONGRESS_GDN_F356W.yaml
"""

from __future__ import annotations

import argparse
import collections
import os
import sys
import time
from multiprocessing import Pool

import numpy as np
import matplotlib
matplotlib.use("Agg")

from astropy.io import fits, ascii
from astropy.table import Table

from nircam_wfss.config import PipelineConfig
from nircam_wfss.background import (
    background_grism_stage2,
    assignwcs_grism_stage2,
    my_grism_bkg_subtraction,
)
from nircam_wfss.preprocessing import grism_hot_pix_rejection, my_grism_cont_subtraction
from nircam_wfss.imaging import reduce_img_stage2, my_daofind_sw_fits, calibrate_astrometry
from nircam_wfss.dispersion import GrismConf
from nircam_wfss.extraction import (
    extract_2d_spec_worker,
    extract_1d_spec_worker,
    extract_2d_emline_worker_drizzle,
    extract_2d_emline_worker_simple,
)
from nircam_wfss.pom import build_POM_applied_catalog

import warnings
warnings.filterwarnings("ignore", message="'obsfix' made the change")
warnings.filterwarnings("ignore", message="'datfix' made the change")
warnings.filterwarnings("ignore", message = "Card is too long")
warnings.filterwarnings("ignore", message = "divide by zero")
warnings.filterwarnings("ignore", message = "invalid value encountered")
warnings.filterwarnings("ignore", message = "Input data contains invalid values")
warnings.filterwarnings("ignore", message = "Input association file contains path information")
warnings.filterwarnings("ignore", message="Double sampling check FAILED")
warnings.filterwarnings("ignore", message="All-NaN slice encountered")
warnings.filterwarnings("ignore", message="Mean of empty slice")
warnings.filterwarnings("ignore", message="Degrees of freedom <= 0")
warnings.filterwarnings("ignore", message="Cannot merge meta key")


# ---------------------------------------------------------------------------
# Helpers for pool-compatible partial function calls
# ---------------------------------------------------------------------------

def _assignwcs_worker(args):
    """Unpack (file, v1p5_dir, overwrite) and call assignwcs_grism_stage2."""
    return assignwcs_grism_stage2(*args)


def _bkg_worker(args):
    """Unpack (file, cali_support_dir, plot_dir) and call my_grism_bkg_subtraction."""
    return my_grism_bkg_subtraction(*args)


def _hotpix_worker(args):
    """Unpack (file, sigma, plot_dir) and call grism_hot_pix_rejection."""
    return grism_hot_pix_rejection(*args)


def _cont_worker(args):
    """Unpack (file, plot_dir) and call my_grism_cont_subtraction."""
    return my_grism_cont_subtraction(*args)

def _stage2_img_worker(args):
    """Unpack (file, output_dir, wing, pid) and call reduce_img_stage2."""
    return reduce_img_stage2(*args)

def _daofind_worker(args):
    """Unpack (file, direct_image_dir, astrometry_dir) and call my_daofind_sw_fits."""
    return my_daofind_sw_fits(*args)

def _extract_2d_spec_worker(args):
    """Unpack (source_id, coord, paths, grism_conf, output_dir) and call extract_2d_spec."""
    return extract_2d_spec_worker(*args)

def _extract_1d_spec_worker(args):
    """Unpack (spec2d_path, extraction_dir, do_boxcar, grism_conf, mosaic params) and call extract_1d_spec_worker."""
    return extract_1d_spec_worker(*args)

def _extract_2d_emline_worker_drizzle(args):
    """Unpack args and call extract_2d_emline_worker_drizzle."""
    return extract_2d_emline_worker_drizzle(*args)

def _extract_2d_emline_worker_simple(args):
    """Unpack args and call extract_2d_emline_worker_simple."""
    return extract_2d_emline_worker_simple(*args)


# ---------------------------------------------------------------------------
# Main pipeline function
# ---------------------------------------------------------------------------

def main(argv: list[str] | None = None) -> int:
    """
    Run the full NIRCam WFSS grism extraction pipeline.

    Parses one positional argument: the path to a YAML configuration file.
    Pipeline stages 2a–8 are executed in sequence using the settings from
    that file (see the module-level docstring for a stage overview).

    Parameters
    ----------
    argv:
        Argument list forwarded to ``argparse``.  ``None`` falls back to
        ``sys.argv[1:]``.

    Returns
    -------
    exit_code:
        0 on success, 1 if the configuration file is not found.
    """

    # =========================================================================
    # Parse command-line arguments
    # =========================================================================
    parser = argparse.ArgumentParser(
        prog="run_pipeline.py",
        description="NIRCam WFSS grism extraction pipeline",
    )
    parser.add_argument(
        "config",
        metavar="CONFIG.yaml",
        help="Path to the YAML pipeline configuration file.",
    )
    args = parser.parse_args(argv)

    if not os.path.isfile(args.config):
        print("ERROR: config file not found: %s" % args.config, file=sys.stderr)
        return 1

    cfg = PipelineConfig.from_yaml(args.config)
    print("Loaded configuration from: %s" % args.config)
    print("  PID=%d  filter=%s  n_procs=%d" % (cfg.pid, cfg.grism_filter, cfg.n_procs))

    # CRDS server URL (override if not set in shell profile)
    os.environ.setdefault("CRDS_SERVER_URL", "https://jwst-crds.stsci.edu")
    print("My HOME directory is  :", os.environ['HOME'])
    print("CRDS server at STScI  :", os.environ['CRDS_SERVER_URL'])
    print("CRDS_PATH on computer :", os.environ['CRDS_PATH'])

    # =========================================================================
    # STAGE 2a – Assign WCS + flat-field: rate.fits → lv1.5 files
    # =========================================================================
    print("\n========== STAGE 2a: Assign WCS + flat-field ==========")
    # select the grism files to process 
    list_rate_this_band = cfg.select_rate_files(filter=cfg.grism_filter)
    # assign WCS and flat-field
    wcs_args = [(f, cfg.v1p5_dir, cfg.overwrite) for f in list_rate_this_band]
    with Pool(min(cfg.n_procs, len(list_rate_this_band))) as pool:
        pool.map(_assignwcs_worker, wcs_args)
    # get the v1.5 files 
    list_v1p5_this_band = cfg.select_lv1p5_files(filter=cfg.grism_filter)

    # =========================================================================
    # STAGE 2b – Build median super-sky backgrounds (if not using pre-computed)
    # =========================================================================
    print("\n========== STAGE 2b: Build super-sky backgrounds ==========")
    # Group lv1.5 files by (filter, module, pupil)
    grp_headers = [(fits.getheader(f), f) for f in list_v1p5_this_band]
    groups: dict[tuple, list] = collections.defaultdict(list)
    for hd, f in grp_headers:
        key = (hd["FILTER"], hd["MODULE"], hd["PUPIL"])
        groups[key].append(f)  # load SCI frame

    for (tmp_flt, tmp_mod, tmp_pup), fnames in groups.items():
        background_grism_stage2(
            fnames, tmp_mod, tmp_pup, tmp_flt,
            cali_support_dir=cfg.cali_support_dir,
            overwrite=False,
        )

    # =========================================================================
    # STAGE 2c – Background subtraction + 1/f noise removal (v1.5 inplace)
    # =========================================================================
    print("\n========== STAGE 2c: Background subtraction ==========")
    bkg_args = [(f, cfg.cali_support_dir, cfg.plot_dir) for f in list_v1p5_this_band]
    with Pool(min(cfg.n_procs, len(list_v1p5_this_band))) as pool:
        pool.map(_bkg_worker, bkg_args)
    # =========================================================================
    # STAGE 2d – Hot-pixel rejection (v1.5 inplace)
    # =========================================================================
    print("\n========== STAGE 2d: Hot-pixel rejection ==========")
    hotpix_args = [(f, cfg.sigma_hot, cfg.plot_dir) for f in list_v1p5_this_band]
    with Pool(min(cfg.n_procs, len(list_v1p5_this_band))) as pool:
        pool.map(_hotpix_worker, hotpix_args)

    # =========================================================================
    # STAGE 2e – Continuum subtraction → EMLINE extensions (v1.5 append)
    # =========================================================================
    print("\n========== STAGE 2e: Continuum subtraction ==========")
    cont_args = [(f, cfg.plot_dir) for f in list_v1p5_this_band]
    with Pool(min(cfg.n_procs, len(list_v1p5_this_band))) as pool:
        pool.map(_cont_worker, cont_args)

    # =========================================================================
    # STAGE 3 – Reduce SW direct images (calwebb_image2)
    # =========================================================================
    print("\n========== STAGE 3: SW direct image Stage-2 reduction ==========")
    if len(cfg.list_rate_sw) > 0:
        sw_args = [
            (f, cfg.direct_image_dir, False)
            for f in cfg.list_rate_sw
        ]
        with Pool(min(cfg.n_procs, len(cfg.list_rate_sw))) as pool:
            pool.map(_stage2_img_worker, sw_args)
        print("%d SW cal files processed." % len(cfg.list_rate_sw))
    else:
        print("No SW rate files found – skipping Stage-3.")

    # =========================================================================
    # STAGE 4 – Astrometry Calibration
    # =========================================================================
    print("\n========== STAGE 4: Astrometry Calibration ==========")
    if os.path.isfile(cfg.astrometry_cal_table):
        tb_sw_astrometry = ascii.read(cfg.astrometry_cal_table)
        print("Loaded astrometry table with %d entries." % len(tb_sw_astrometry))
    else:
        # calibrate astrometry by runing DAOFIND & cross-match to reference catalog.
        daofind_args = [
            (f, cfg.direct_image_dir, cfg.astrom_dir) for f in cfg.list_rate_sw
        ]
        with Pool(min(cfg.n_procs, len(cfg.list_rate_sw))) as pool:
            pool.map(_daofind_worker, daofind_args)
        print("finished Pool -- beginning calibrate astrometry")
        tb_sw_astrometry = calibrate_astrometry(
            cfg.list_rate_sw,
            cfg.direct_image_dir,
            cfg.v1p5_dir,
            cfg.astrom_dir,
            cfg.astrometry_cal_table,
            cfg.astrometry_ref_table,
            False,
        )

    # =========================================================================
    # STAGE 5 – Load source catalogue & build POM-applied catalog per frame
    # =========================================================================
    print("\n========== STAGE 5: Prepare POM-applied catalog ==========")
    # Load source catalog
    if not os.path.isfile(cfg.source_catalog_path):
        raise FileNotFoundError(
            "Source catalogue %s not found." % cfg.source_catalog_path
        )
    tb_source = Table.read(cfg.source_catalog_path)
    print("Loaded %d sources from %s." % (len(tb_source), os.path.basename(cfg.source_catalog_path)))
    # Build POM-applied catalog for each grism frame, and group by source ID
    frame_ID, frame_path = build_POM_applied_catalog(
        tb_source, tb_sw_astrometry, 
        list_v1p5_this_band,
        cali_support_dir=cfg.cali_support_dir,
        default_POM_trans_dir=cfg.default_POM_trans_dir,
        POM_catalog_dir=cfg.POM_catalog_dir,
        dx=0.0, dy=0.0, POM_threshold=0.0
    )

    # =========================================================================
    # STAGE 6 – Extract 2D spectra for each source
    # =========================================================================
    print("\n========== STAGE 6: Extract 2D Grism Spectra ==========")
    # Load grism configuration
    grism_conf = GrismConf(filter=cfg.grism_filter, config=cfg)
    # Extract 2D spectra for each source in the catalog. 
    spec2d_args = []
    for i in range(len(tb_source)):
        objid = str(tb_source["ID"][i])
        POM_of_this_source = frame_path["path_"+objid]
        if len(POM_of_this_source) == 0:
            print("Warning: no POM catalog found for source ID %s in any grism frame. Skipping." % objid)
            continue
        spec2d_args.append((
            frame_ID[objid],
            frame_path["path_"+objid],
            list_v1p5_this_band,
            tb_source[i],
            grism_conf,
            cfg.aperture_pix,
            cfg.grism_filter,
            cfg.extract_dir,
            cfg.bunit_spec2d,
            True,
            cfg.overwrite_spec2d,
        ))
    with Pool(min(cfg.n_procs, len(tb_source))) as pool:
        pool.map(_extract_2d_spec_worker, spec2d_args)


    # =========================================================================
    # STAGE 7 – Extract 1D spectra for each source
    # =========================================================================
    print("\n========== STAGE 7: Extract 1D Grism Spectra ==========")
    spec1d_args = []
    for i in range(len(tb_source)):
        spec2d_path = os.path.join(
            cfg.extract_dir,
            "spec_2d_%s_ID%s_allcoadd.fits" % (grism_conf.filter, tb_source["ID"][i])
        )
        if not os.path.isfile(spec2d_path):
            continue
        spec1d_args.append((
            spec2d_path,
            cfg.extract_dir,
            True,
            grism_conf,
            cfg.image_mosaic_dir,
            cfg.image_mosaic_filename_fmt,
            cfg.image_mosaic_rgb_bands,
            cfg.image_mosaic_field,
            cfg.plot_dir,
        ))
    with Pool(min(cfg.n_procs, len(tb_source))) as pool:
        pool.map(_extract_1d_spec_worker, spec1d_args)

    # =========================================================================
    # STAGE 8 – Extract 2D emission-line cutouts for each source
    # =========================================================================
    # Requires 'z_grism' and 'name_line_exp' columns in the source catalogue.
    print("\n========== STAGE 8: Extract 2D Emission-Line Cutouts ==========")
    if "z_grism" in tb_source.colnames and "name_line_exp" in tb_source.colnames:
        if cfg.coadd_method == "drizzle":
            print("Using drizzle-based extraction for emission-line cutouts.")
            emline_args = []
            for i in range(len(tb_source)):
                objid = str(tb_source["ID"][i])
                POM_of_this_source = frame_path["path_" + objid]
                if len(POM_of_this_source) == 0:
                    continue
                emline_args.append((
                    frame_ID[objid],
                    frame_path["path_" + objid],
                    list_v1p5_this_band,
                    tb_source[i],
                    grism_conf,
                    cfg.grism_filter,
                    cfg.extract_dir,
                    cfg.cutout_size_drizzle,      # cutout_size
                    cfg.finalscale_drizzle,
                    cfg.pixfrac_drizzle,
                    cfg.psf_oversample,  # PSF super-sampling factor
                ))
            if emline_args:
                with Pool(min(cfg.n_procs, len(emline_args))) as pool:
                    pool.map(_extract_2d_emline_worker_drizzle, emline_args)
        else:
            print("Using simple direct cutout extraction for emission-line cutouts.")
            emline_args = []
            for i in range(len(tb_source)):
                objid = str(tb_source["ID"][i])
                POM_of_this_source = frame_path["path_" + objid]
                if len(POM_of_this_source) == 0:
                    continue
                emline_args.append((
                    tb_source[i],
                    grism_conf,
                    cfg.grism_filter,
                    cfg.extract_dir,
                    cfg.cutout_size_simple,
                    cfg.psf_oversample,
                ))
            if emline_args:
                with Pool(min(cfg.n_procs, len(emline_args))) as pool:
                    pool.map(_extract_2d_emline_worker_simple, emline_args)
    else:
        print("  Skipping: source catalogue lacks 'z_grism' or 'name_line_exp' columns.")

    return 0
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    sys.exit(main())
