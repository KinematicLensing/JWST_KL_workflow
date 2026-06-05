"""
nircam_wfss
===========
NIRCam Wide-Field Slitless Spectroscopy (WFSS) grism extraction pipeline.

Converted from the pedagogical Jupyter notebook
  NIRCam_grism_extraction_code_example_FSun.ipynb
  by Fengwu Sun (CfA/Harvard).

Package layout
--------------
config.py       PipelineConfig dataclass; emission-line wavelength table; constants.
noise.py        1/f noise subtraction; robust background estimation.
wcs_utils.py    Attitude matrix; astrometric correction helpers.
dispersion.py   Grism dispersion/trace polynomials; grism_conf_preparation.
pom.py          POMData class; pick-off mirror vignetting checks.
background.py   Median super-sky creation; WCS+flat-field (lv1.5); bkg subtraction.
preprocessing.py Hot-pixel rejection; continuum/emission-line separation.
imaging.py      Stage-2 direct-image reduction; SW DAOStarFinder astrometry.
extraction.py   2-D spectral extraction; multi-frame co-addition; resampling.
plotting.py     Matplotlib corner-label utilities.
run_pipeline.py Stand-alone pipeline script.
"""

from .config import PipelineConfig, EML_LAB, NIRCAM_LW_ZP, COSMO

from .noise import my_1overf_subtraction, robust_median_bkg

from .wcs_utils import rotate, attitude, apply_astrometric_correction

from .dispersion import (
    fit_disp_order32,
    fit_disp_order23,
    grism_conf_preparation,
    linear,
    gauss,
    gauss_cont_prof,
)

from .pom import POMData, is_pickoff, is_pickoff_ps

from .background import (
    get_crds_dict_from_fits_header,
    background_grism_stage2,
    assignwcs_grism_stage2,
    my_grism_bkg_subtraction,
)

from .preprocessing import grism_hot_pix_rejection, my_grism_cont_subtraction

from .imaging import reduce_img_stage2, my_daofind_sw_fits

from .extraction import extract_2d_spec, store_all_2d_spec, resample_spec2d_wmin_wmax

from .plotting import get_corner_pos, corner_text

__all__ = [
    # config
    "PipelineConfig", "EML_LAB", "NIRCAM_LW_ZP", "COSMO",
    # noise
    "my_1overf_subtraction", "robust_median_bkg",
    # wcs
    "rotate", "attitude", "apply_astrometric_correction",
    # dispersion
    "fit_disp_order32", "fit_disp_order23", "grism_conf_preparation",
    "linear", "gauss", "gauss_cont_prof",
    # pom
    "POMData", "is_pickoff", "is_pickoff_ps",
    # background
    "get_crds_dict_from_fits_header",
    "background_grism_stage2", "assignwcs_grism_stage2", "my_grism_bkg_subtraction",
    # preprocessing
    "grism_hot_pix_rejection", "my_grism_cont_subtraction",
    # imaging
    "reduce_img_stage2", "my_daofind_sw_fits",
    # extraction
    "extract_2d_spec", "store_all_2d_spec", "resample_spec2d_wmin_wmax",
    # plotting
    "get_corner_pos", "corner_text",
]
