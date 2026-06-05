"""
pom.py
------
NIRCam Pick-Off Mirror (POM) vignetting utilities.

The POM is a pick-off mirror that partially blocks the field of view for
NIRCam LW grism observations.  Sources that fall behind the POM will have
their spectra partially or completely vignetted.

Two levels of information are provided:

1. ``POMData`` – loads and caches the POM transmission maps for both modules.
   Pass one instance to functions instead of module-level globals.

2. ``is_pickoff`` / ``is_pickoff_ps`` – check whether a source at a given
   pixel position is obscured by the POM or yields only a partial spectrum.
"""

from __future__ import annotations

import os

import numpy as np
from astropy.io import fits
from astropy import wcs
from astropy.table import Table
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import ascii


class POMData:
    """
    Container for NIRCam LW POM transmission maps.

    Loads the FITS files once on construction and exposes the data arrays
    and pixel offsets for both Module A and Module B.

    Parameters
    ----------
    pom_data_dir:
        Directory containing ``NIRCAM_LW_POM_ModA_trans.fits`` and
        ``NIRCAM_LW_POM_ModB_trans.fits``.

    Notes
    -----
    The FITS files can be downloaded from:
    https://github.com/npirzkal/GRISM_NIRCAM/tree/master/V4
    """

    def __init__(self, pom_data_dir: str = "./data/GRISM_NIRCAM") -> None:
        self.pom_data_dir = pom_data_dir

        path_a = os.path.join(pom_data_dir, "NIRCAM_LW_POM_ModA_trans.fits")
        path_b = os.path.join(pom_data_dir, "NIRCAM_LW_POM_ModB_trans.fits")

        with fits.open(path_a) as hdul_a:
            self.trans_A: np.ndarray = hdul_a[1].data
            self.xy_start_A: tuple[int, int] = (
                int(hdul_a[1].header["NOMXSTRT"]),
                int(hdul_a[1].header["NOMYSTRT"]),
            )

        with fits.open(path_b) as hdul_b:
            self.trans_B: np.ndarray = hdul_b[1].data
            # Note: uses ModA Y-start intentionally (matches original notebook)
            self.xy_start_B: tuple[int, int] = (
                int(hdul_b[1].header["NOMXSTRT"]),
                int(hdul_a[1].header["NOMYSTRT"]),  # intentional – see notebook
            )

    def get(self, module: str) -> tuple[tuple[int, int], np.ndarray]:
        """Return ``(xy_start, trans_array)`` for the given module ('A' or 'B')."""
        module = module.upper()
        if module == "A":
            return self.xy_start_A, self.trans_A
        elif module == "B":
            return self.xy_start_B, self.trans_B
        raise ValueError("module must be 'A' or 'B'")


# ---------------------------------------------------------------------------
# POM vignetting checks
# ---------------------------------------------------------------------------

def is_pickoff(
    x: int | float | list | np.ndarray,
    y: int | float | list | np.ndarray,
    module: str = "A",
    pom_data: POMData | None = None,
    pom_data_dir: str = "./data/GRISM_NIRCAM",
) -> np.ndarray:
    """
    Determine whether source(s) at pixel (x, y) are obscured by the pick-off mirror.

    Parameters
    ----------
    x, y:
        Pixel coordinate(s) in the NIRCam detector frame.
    module:
        Detector module: ``'A'`` or ``'B'``.
    pom_data:
        Pre-loaded ``POMData`` instance.  If ``None``, one is constructed
        on the fly from ``pom_data_dir``.
    pom_data_dir:
        Used only when ``pom_data`` is ``None``.

    Returns
    -------
    arr_trans:
        Integer array (same length as inputs).  0 = obscured, 1 = clear.
    """
    if pom_data is None:
        pom_data = POMData(pom_data_dir)

    xy_start, pom_trans = pom_data.get(module)

    tmp_x = np.atleast_1d(np.int32(x))
    tmp_y = np.atleast_1d(np.int32(y))
    if len(tmp_x) != len(tmp_y):
        raise ValueError("x and y must have the same length")

    tmp_x = tmp_x + xy_start[0]
    tmp_y = tmp_y + xy_start[1]

    arr_trans = []
    ny, nx = pom_trans.shape
    for j in range(len(tmp_x)):
        if tmp_x[j] < 0 or tmp_y[j] < 0 or tmp_x[j] > nx - 1 or tmp_y[j] > ny - 1:
            arr_trans.append(0)
        else:
            arr_trans.append(int(pom_trans[tmp_x[j], tmp_x[j]]))

    return np.array(arr_trans, dtype=int)


def is_pickoff_ps(
    x: int | float | list | np.ndarray,
    y: int | float | list | np.ndarray,
    module: str = "A",
    grism_filter: str = "F444W",
    pupil: str = "R",
    spec_cov_data_dir: str = "./data/FSun_cal",
) -> np.ndarray:
    """
    Check whether source(s) yield complete, partial, or no spectra due to the POM. PS = partial spectrum

    Uses the spectral-coverage FITS files (``FSun_SpecCov_<filter>_<module>_<pupil>.fits``)
    to account for filter-dependent spectral completeness.

    Parameters
    ----------
    x, y:
        Pixel coordinate(s) in the NIRCam detector frame.
    module:
        Detector module: ``'A'`` or ``'B'``.
    grism_filter:
        NIRCam filter (e.g. ``'F444W'``, ``'F356W'``).
    pupil:
        Grism pupil: ``'R'`` (GRISMR) or ``'C'`` (GRISMC).
    spec_cov_data_dir:
        Directory containing the spectral-coverage FITS files.
        Files are available at https://magnif.as.arizona.edu/~fsun/home/data/spec_cov_fits/

    Returns
    -------
    arr_trans:
        Integer array with spectral coverage fraction (0–100 %).
        0 = completely obscured, 100 = full spectrum.
    """
    path_ps = os.path.join(
        spec_cov_data_dir,
        "FSun_SpecCov_%s_%s_%s.fits" % (grism_filter, module, pupil),
    )
    if not os.path.isfile(path_ps):
        raise FileNotFoundError("Spectral coverage file not found: %s" % path_ps)

    with fits.open(path_ps) as hdul:
        xy_start = (
            int(hdul[0].header["NOMXSTRT"]),
            int(hdul[0].header["NOMYSTRT"]),
        )
        pom_trans = hdul[0].data

    tmp_x = np.atleast_1d(np.int32(x))
    tmp_y = np.atleast_1d(np.int32(y))
    if len(tmp_x) != len(tmp_y):
        raise ValueError("x and y must have the same length")

    tmp_x = tmp_x + xy_start[0]
    tmp_y = tmp_y + xy_start[1]

    arr_trans = []
    ny, nx = pom_trans.shape
    for j in range(len(tmp_x)):
        if tmp_x[j] < 0 or tmp_y[j] < 0 or tmp_x[j] > nx - 1 or tmp_y[j] > ny - 1:
            arr_trans.append(0)
        else:
            arr_trans.append(int(pom_trans[tmp_y[j], tmp_x[j]]))

    return np.array(arr_trans, dtype=int)


def build_POM_applied_catalog(
    tb_source: "astropy.table.Table",  
    tb_sw_astrometry: "astropy.table.Table",
    all_lv1p5_list: list[str],
    cali_support_dir: str = "./data/FSun_cal",
    default_POM_trans_dir: str = "./data/GRISM_NIRCAM",
    POM_catalog_dir: str = "./grism_demo/POM_catalog",
    dx=0.0,
    dy=0.0,
    POM_threshold: float = 0.05,
):
    """
    Build a copy of the source catalogue with POM spectral coverage applied.

    Adds a new column 'POM_SPEC_COV' with values from 0 to 100 indicating the
    percentage of the spectrum that is expected to be observable given the
    source's pixel position and the POM vignetting pattern.

    Parameters
    ----------
    tb_source:
        Source catalogue with columns 'x' and 'y' for pixel positions. Note that
        for different grism filters, the source catalog should be different. 
    tb_sw_astrometry:
        Table containing astrometric offsets measured from SW direct images.
    all_lv1p5_list:
        List of paths to the level-1.5 grism exposure FITS files.
    cali_support_dir:
        Directory containing the spectral-coverage FITS files for partial spectrum checks.
    default_POM_trans_dir:
        Directory containing the POM coverage data
    dx, dy:
        Pixel offsets to apply to the source positions due to filter offset. 
    POM_catalog_dir:
        Path to save the POM-applied catalog.
    POM_threshold:
        Minimum spectral coverage fraction (0–100) for a source to be considered observable.
        If set to zero, the function will call the POM coverage rather than POM spectral coverage fraction.
    Returns
    -------
    grism_frame_ID_per_source:
        Dictionary mapping source ID to list of grism frame indices where the source is observable.
    grism_frame_path_per_source:
        Dictionary mapping source ID to list of POM-applied catalog paths where the source is observable.
    """
    POM_catalogs = []
    POM_catalog_paths = []
    for i, tmp_lv1p5_path in enumerate(all_lv1p5_list):
        # Apply astrometric offsets measured from simultaneous SW direct images.
        # Pointing accuracy can vary; the SW-derived dRA/dDEC/theta correction
        # is applied to each grism frame's SCI header before projecting sources.
        tmp_rate_path_base = tmp_lv1p5_path.split('/')[-1].split('long_rate')[0]
        print('[%3d]' % i, tmp_rate_path_base)

        # Load header for metadata
        tmp_grism_hd_1st = fits.getheader(tmp_lv1p5_path, 0)
        tmp_grism_hd_sci = fits.getheader(tmp_lv1p5_path, 'sci')
        tmp_filter = tmp_grism_hd_1st['filter']
        
        source_coords = SkyCoord(tb_source["RA"], tb_source["DEC"], unit=(u.deg, u.deg))

        item_sw_astrom = tb_sw_astrometry[tb_sw_astrometry['expName'] == tmp_rate_path_base]
        if len(item_sw_astrom) > 0:
            ## if there is astrometric information, consider that
            item_sw_astrom = item_sw_astrom[0]
            tmp_grism_hd_sci['CRVAL1'] -= item_sw_astrom['dRA'] / np.cos(np.deg2rad(tmp_grism_hd_sci['CRVAL2'])) / 3600.
            tmp_grism_hd_sci['CRVAL2'] -= item_sw_astrom['dDEC'] / 3600. 
            tmp_cd_matrix = np.array([[tmp_grism_hd_sci['CD1_1'], tmp_grism_hd_sci['CD1_2']],
                                    [tmp_grism_hd_sci['CD2_1'], tmp_grism_hd_sci['CD2_2']]])
            tmp_rot_matrix = np.array([[np.cos(np.deg2rad(-item_sw_astrom['theta'])), - np.sin(np.deg2rad(-item_sw_astrom['theta']))],
                                    [np.sin(np.deg2rad(-item_sw_astrom['theta'])),   np.cos(np.deg2rad(-item_sw_astrom['theta']))]])
            tmp_cd_matrix = np.matmul(tmp_cd_matrix, tmp_rot_matrix)
            tmp_grism_hd_sci['CD1_1'] = tmp_cd_matrix[0][0]
            tmp_grism_hd_sci['CD1_2'] = tmp_cd_matrix[0][1]
            tmp_grism_hd_sci['CD2_1'] = tmp_cd_matrix[1][0]
            tmp_grism_hd_sci['CD2_2'] = tmp_cd_matrix[1][1]
        else:
            ## if no astrometric information avaialble, skip this step and directly use RA and DEC
            print('Warning: No SW images found for %s, skip astrom correction! ' % tmp_rate_path_base)
        
        # Project source RA/DEC into grism-frame pixel coordinates.
        tmp_grism_wcs = wcs.WCS(tmp_grism_hd_sci)
        idx_this_field = np.where(
            (np.abs((tb_source["RA"] - tmp_grism_hd_sci['crval1']) * np.cos(np.deg2rad(tmp_grism_hd_sci['crval2']))) < 4 / 60.) &
            (np.abs(tb_source["DEC"] - tmp_grism_hd_sci['crval2']) < 4 / 60.))[0]
        if len(idx_this_field) == 0: continue
        pixelx, pixely  = wcs.utils.skycoord_to_pixel(source_coords[idx_this_field], tmp_grism_wcs)
        # apply filter offset in needed
        pixelx, pixely = pixelx + dx, pixely + dy

        tb_sub = tb_source[idx_this_field]

        if 'F444W_mag' in tb_sub.colnames: tmp_mag_auto = tb_sub['F444W_mag'].data
        elif 'F200W_mag' in tb_sub.colnames: tmp_mag_auto = tb_sub['F200W_mag'].data
        tb_pom_applied = Table(data = [tb_sub['ID'].data, tb_sub["RA"], tb_sub["DEC"], 
                                    pixelx, pixely, tmp_mag_auto,],
                    names = ['Index', 'ra', 'dec', 'pixel_x', 'pixel_y', 'MAG_AUTO'])
        # whether the source falls on POM or not
        # If POM_threshold is set to zero, we only care if the source is picked by the POM, we don't care how much spectrum it yield
        # If POM_threshold is set to a value between 0 and 1, we require the spectrum to be at least POM_threshold complete
        # People can set POM_threshold to zero first, and then reject sources don't yield meaningful spectra. 
        if POM_threshold > 1e-5:
            arg_is_pickoff = np.where(is_pickoff_ps(pixelx, pixely, 
                    module = tmp_grism_hd_1st['module'], 
                    grism_filter = tmp_grism_hd_1st['filter'],
                    pupil = tmp_grism_hd_1st['pupil'][-1],
                    spec_cov_data_dir = cali_support_dir) > POM_threshold)[0]
        else:
            arg_is_pickoff = np.where(is_pickoff(pixelx, pixely, 
                    module = tmp_grism_hd_1st['module'],
                    pom_data_dir = default_POM_trans_dir))[0]
        tb_pom_applied = tb_pom_applied[arg_is_pickoff]
        tb_pom_applied['ra'].info.format = '.6f'
        tb_pom_applied['dec'].info.format = '.6f'
        tb_pom_applied['pixel_x'].info.format = '.3f'
        tb_pom_applied['pixel_y'].info.format = '.3f'
        tb_pom_applied['MAG_AUTO'].info.format = '.3f'
        try:
            tb_pom_applied['MAGERR_AUTO'].info.format = '.3f'
            tb_pom_applied['A_IMAGE'].info.format = '.3f'
            tb_pom_applied['B_IMAGE'].info.format = '.3f'
            tb_pom_applied['THETA_IMAGE'].info.format = '.2f'
            tb_pom_applied['CLASS_STAR'].info.format = '.3f'
        except KeyError: pass
        
        path_tb_pom_applied = os.path.join(POM_catalog_dir, 
            tmp_rate_path_base + '_%s_dirimg_sources.list' % tmp_filter)
        if i == 0:
            if os.path.isdir(os.path.dirname(path_tb_pom_applied)) == False:
                os.mkdir(os.path.dirname(path_tb_pom_applied))
        ascii.write(tb_pom_applied, path_tb_pom_applied, overwrite = True)
        POM_catalogs.append(tb_pom_applied)
        POM_catalog_paths.append(path_tb_pom_applied)

    # Group POM-applied catalog by source ID for later per-source extraction.
    # ID of exposure for each source
    grism_frame_ID_per_source = dict()
    # path of exposure (*/jw*.fits) for each source
    grism_frame_path_per_source = dict()
    
    for idx in tb_source["ID"][:]:
        grism_frame_ID_per_source['%s' % idx] = []
        grism_frame_path_per_source['path_%s' % idx] = []
    # For each object in each POM-applied catalog, register 
    # 1. the exposure index in the file list
    # 2. the path of the POM-applied catalog 
    for i, tmp_sl_table in enumerate(POM_catalogs):
        for idx in tmp_sl_table['Index'].data:
            grism_frame_ID_per_source['%s' % idx].append(i)
            grism_frame_path_per_source['path_%s' % idx].append(POM_catalog_paths[i])
    N_frame_per_source = np.array([len(grism_frame_ID_per_source[x]) \
                                   for x in grism_frame_ID_per_source.keys()])
    print('%d sources may yield spectra' % (np.sum(N_frame_per_source > 0)))

    return grism_frame_ID_per_source, grism_frame_path_per_source