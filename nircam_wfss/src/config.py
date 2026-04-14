"""
config.py
---------
Pipeline configuration dataclass and physical constants shared across all modules.

Design notes
~~~~~~~~~~~~
All directory paths and run-time tuning parameters live in ``PipelineConfig``.
Pass a ``PipelineConfig`` instance into every function that previously relied on
module-level globals (``calibrated_dir``, ``extract_dir``, ``n_procs``, etc.).
This makes every function independently testable and avoids hidden state.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from typing import Optional, Any

from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
from astropy.io import fits
import numpy as np
import glob as _glob

# ---------------------------------------------------------------------------
# Physical / observational constants
# ---------------------------------------------------------------------------

#: Default NIRCam LW zeropoint  (AB mag, -2.5 log10 of flux conversion factor)
NIRCAM_LW_PIXSCALE: float = 0.0629  # arcsec/pixel
NIRCAM_LW_ZP: float = -2.5 * np.log10((u.MJy / u.sr * (NIRCAM_LW_PIXSCALE*u.arcsec)**2 / (3631 * u.Jy)).cgs.value)

#: Default cosmology used throughout the pipeline
COSMO = FlatLambdaCDM(H0=70, Om0=0.3)

#: Emission line rest-frame wavelengths in microns.
#: Keys are common spectroscopic identifiers.
EML_LAB: dict[str, float] = {
    # Hydrogen Brackett series
    "BrA":  4.051, "BrB":  2.625, "BrG":  2.166, "BrD":  1.944,
    # Hydrogen Paschen series
    "PaA":  1.875, "PaB":  1.282, "PaG":  1.094,
    # Helium
    "HeI":  1.083,
    # Hydrogen Balmer series
    "Ha":   0.6563, "Hb":  0.4861, "Hg":  0.4341, "Hd":  0.4102,
    # Oxygen forbidden lines
    "[OII]3727":  0.3727, "[OII]3729":  0.3729,
    "[OIII]4959": 0.4959, "[OIII]5007": 0.5007,
    "[OI]6300":   0.6300,
    # Nitrogen forbidden lines
    "[NII]6548":  0.6548, "[NII]6583":  0.6583,
    # Sulfur forbidden lines
    "[SII]6717":  0.6717, "[SII]6731":  0.6731,
    # Carbon / Silicon / Magnesium UV lines
    "CIV1548":    0.1548, "CIV1551":    0.1551,
    "CIII]1907":  0.1907, "CIII]1909":  0.1909,
    "MgII2796":   0.2796, "MgII2803":   0.2803,
    # Lyman series
    "Lya":  0.1216,
    # Molecular hydrogen
    "H2_1-0S1": 2.1218,
    # PAH features
    "PAH3.3":   3.300,
    # CO bandheads
    "CO(2-0)":  2.294,
    # Mid-IR lines
    "NeII12.8": 12.81, "NeIII15.6": 15.56,
    "ArII6.99": 6.985, "ArIII9.0":  8.991,
    # Ionic forbidden lines
    "[FeII]1.257": 1.257, "[FeII]1.644": 1.644,
}


# ---------------------------------------------------------------------------
# Pipeline configuration
# ---------------------------------------------------------------------------

@dataclass
class PipelineConfig:
    """
    Central configuration object for the NIRCam WFSS pipeline.

    Pass one instance of this class to every pipeline function.  Avoids the
    use of module-level global variables for paths and shared settings.

    Parameters
    ----------
    pid: 
        Programme ID (cosmetic; used in Stage-2 imaging association files).
    data_dir:
        Root directory that contains all raw and calibrated data sub-folders.
    calibrated_dir:
        Directory where level-1.5 (WCS-assigned, flat-fielded) grism FITS files
        and median super-sky backgrounds are stored.  Created automatically if
        it does not exist.
    extract_dir:
        Output directory for extracted 2-D spectra.  Created automatically if
        it does not exist.
    direct_image_dir:
        Directory containing NIRCam SW/LW direct-imaging *rate* and *cal* files.
    grism_data_dir:
        Directory containing the original NIRCam LW grism *rate* files (stage-1
        outputs from the JWST pipeline).
    cali_support_dir:
        Directory with calibration supportive data files, including 
        - spectral-coverage (partial-POM) FITS files produced by F. Sun 
        (``FSun_SpecCov_<filter>_<module>_<pupil>.fits``).
        - grism sensitivity files
        - grism dispersion and spectral tracing solutions
        - median super-sky background images for each filter/pupil/module combination
        - astrometry calibration tables for SW-based astrometric correction
    default_POM_trans_dir:
        Directory to look for the default POM transformation files if not found in 
        cali_support_dir.
        These files are used to determine which sources fall on the POM and should be
          excluded from
    grism_filter:
        NIRCam filter used for the grism observation (e.g. ``'F444W'``).
    n_procs:
        Number of parallel worker processes for pool-based steps.
    # extract_mode:
    #     Which grism pupils to include: ``'all'``, ``'R'``, or ``'C'``.
    # mod_mode:
    #     Which detector modules to include: ``'comb'`` (A+B), ``'modA'``, or
    #     ``'modB'``.
    image_mosaic_dir:
        Directory containing large JWST image mosaics on disk (used for
        direct-image cutouts instead of querying MAST).
    image_mosaic_filename_fmt:
        ``str % (field, band)`` format string that resolves to a mosaic
        FITS filename inside *image_mosaic_dir*.  The first placeholder is
        the field name (e.g. ``'goods-s'``) and the second is the filter
        name in lower-case (e.g. ``'f444w'``).
    image_mosaic_rgb_bands:
        Three-element list ``[blue_band, green_band, red_band]`` giving the
        NIRCam filter names for the blue, green, and red (long-wavelength)
        channels of the diagnostic RGB direct-image thumbnail.
        Default: ``['F090W', 'F200W', 'F444W']``.
    image_mosaic_field:
        Field identifier string inserted into *image_mosaic_filename_fmt*
        (e.g. ``'goods-s'`` or ``'goods-n'``).
    aperture_pix:
        Cross-dispersion aperture half-width in pixels used during 2-D
        spectral extraction.
    sigma_hot:
        Sigma threshold for hot-pixel detection in stage-2d.  Higher values are
        more conservative (fewer hot pixels flagged).  Tune this parameter based
        on the noise properties of your data and the desired balance between
        hot-pixel rejection and preservation of real sources.
    """

    # Directories
    source_catalog_path: str = None
    data_dir:          str = "."
    calibrated_dir:    str = field(default="")
    extract_dir:       str = field(default="")
    direct_image_dir:  str = field(default="")
    direct_image_filename_fmt: str = "jw*_nrc[a-b][1-4]_rate.fits"
    grism_data_dir:    str = field(default="")
    grism_filename_fmt: str = "jw*nrc[ab]long_rate.fits"
    cali_support_dir: str = "./data/FSun_cal"
    default_POM_trans_dir: str = "./data/GRISM_NIRCAM"
    astrometry_cal_table: str = "./data/FSun_cal/PID_3577_SW_wavecal.dat"
    astrometry_ref_table: str = "./data/FSun_cal/goods_charge_F160W_daofind.cat"

    # Image mosaic parameters (for direct-image cutouts read from local disk)
    image_mosaic_dir:          str  = ""
    image_mosaic_filename_fmt: str  = "hlsp_jades_jwst_nircam_%s_%s_v5.0_drz.fits"
    image_mosaic_rgb_bands:    list = field(default_factory=lambda: ["F090W", "F200W", "F444W"])
    image_mosaic_field:        str  = "goods-s"

    # Observation parameters
    pid:           int   = 1895
    grism_filter:  str   = "F444W"
    n_procs:       int   = 4
    # extract_mode:  str   = "all"   # 'all' | 'R' | 'C'
    # mod_mode:      str   = "comb"  # 'comb' | 'modA' | 'modB'
    aperture_pix:  float = 15.0
    sigma_hot:     float = 20.0
    overwrite:     bool = False
    psf_oversample: int  = 4

    def __post_init__(self) -> None:
        # Set default sub-directory names relative to data_dir when not given
        if not self.calibrated_dir:
            self.calibrated_dir = os.path.join(self.data_dir, "grism_cal")
        if not self.extract_dir:
            self.extract_dir = os.path.join(self.data_dir, "extract_2d")
        if not self.direct_image_dir:
            self.direct_image_dir = os.path.join(self.data_dir, "direct_imaging")
        if not self.grism_data_dir:
            self.grism_data_dir = os.path.join(self.data_dir, "grism_raw")
        self.plot_dir = os.path.join(self.calibrated_dir, "plots")
        self.astrom_dir = os.path.join(self.calibrated_dir, "astrom")
        self.POM_catalog_dir = os.path.join(self.calibrated_dir, "POM_catalog")
        self.v1p5_dir = os.path.join(self.calibrated_dir, "lv1p5")
        self.make_dirs()

        # List all rate files in the grism data directory for later use
        self.list_rate_all = np.array(
            sorted(_glob.glob(os.path.join(self.grism_data_dir, 
                                           self.grism_filename_fmt)))
        )
        if len(self.list_rate_all) == 0:
            raise FileNotFoundError(
                "No grism rate files found in: %s" % self.grism_data_dir
            )
        
        # get the filter, pupil, and module for each grism exposure rate file
        self.list_filter,  self.list_pupil, self.list_module = [], [], []
        self.list_target, self.list_exptime = [], []
        for path in self.list_rate_all:
            hd = fits.getheader(path)
            self.list_filter.append(hd["FILTER"])
            self.list_pupil.append(hd["PUPIL"])
            self.list_module.append(hd["MODULE"])
            self.list_target.append(hd["TARGPROP"])
            self.list_exptime.append(hd["EFFEXPTM"])
        self.list_filter = np.array(self.list_filter)
        self.list_pupil  = np.array(self.list_pupil)
        self.list_module = np.array(self.list_module)
        self.list_target = np.array(self.list_target)
        self.list_exptime = np.array(self.list_exptime)
        self.list_v1p5_all = np.array([
                os.path.join(
                    self.v1p5_dir,
                    os.path.basename(f).replace("rate.fits", "rate_lv1.5.fits"),
                )
                for f in self.list_rate_all
            ])
        
        # list all the direct imaging files in SW channel
        self.list_rate_sw = np.array(
            sorted(_glob.glob(
                os.path.join(self.direct_image_dir, 
                self.direct_image_filename_fmt))) if self.direct_image_dir else []
        )

    @classmethod
    def from_yaml(cls, yaml_path: str) -> "PipelineConfig":
        """
        Construct a ``PipelineConfig`` from a YAML configuration file.

        Only keys that correspond to declared dataclass fields are passed to
        the constructor; unrecognised keys are silently ignored so that
        comments or documentation-only entries in the YAML are harmless.

        Parameters
        ----------
        yaml_path:
            Path to the YAML file.  The file is read with
            ``yaml.safe_load``, so no executable Python is ever evaluated.

        Returns
        -------
        PipelineConfig
            Fully initialised configuration object.

        Example
        -------
        .. code-block:: yaml

            pid: 3577
            grism_filter: F356W
            data_dir: /xdisk/timeifler/jiachuanxu/jwst/congress
            n_procs: 20

        .. code-block:: python

            cfg = PipelineConfig.from_yaml("configs/my_run.yaml")
        """
        import yaml
        import dataclasses

        with open(yaml_path) as fh:
            raw: dict[str, Any] = yaml.safe_load(fh) or {}

        # Only forward keys that match declared fields
        field_names = {f.name for f in dataclasses.fields(cls)}
        kwargs = {k: v for k, v in raw.items() if k in field_names}
        return cls(**kwargs)

    def make_dirs(self) -> None:
        """Create output directories if they do not already exist."""
        for d in [self.calibrated_dir,
                  self.plot_dir,
                  self.astrom_dir,
                  self.POM_catalog_dir,
                  self.v1p5_dir,
                  self.extract_dir]:
            os.makedirs(d, exist_ok=True)

    def select_rate_files(self, filter=None) -> None:
        """Select grism rate files matching the target filter and pupil criteria."""
        if filter is None:
            filter = self.grism_filter
        sel_band = (self.list_filter == filter) & (self.list_pupil != "CLEAR")
        list_rate_this_band = self.list_rate_all[sel_band]
        list_rate_this_band.sort()
        print("%d %s grism rate files found." % (len(list_rate_this_band), filter))
        return list_rate_this_band
    
    def select_lv1p5_files(self, filter=None) -> None:
        """Select calibrated lv1.5 files matching the target filter and pupil criteria."""
        if filter is None:
            filter = self.grism_filter
        sel_band = (self.list_filter == filter) & (self.list_pupil != "CLEAR")
        list_lv15_this_band = self.list_v1p5_all[sel_band]
        list_lv15_this_band.sort()
        for f in list_lv15_this_band:
            if not os.path.isfile(f):
                raise FileNotFoundError("Expected lv1.5 file not found: %s" % f)
        print("%d %s grism lv1.5 files found." % (len(list_lv15_this_band), filter))
        return list_lv15_this_band
