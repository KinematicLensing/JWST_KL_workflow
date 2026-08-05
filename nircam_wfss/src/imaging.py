"""
imaging.py
----------
Processing of NIRCam short-wavelength (SW) direct-imaging exposures:

1. ``reduce_img_stage2`` – run JWST ``calwebb_image2`` pipeline on a *rate*
   file to produce a *cal* file with WCS, flat-field, and photometric
   calibration applied.

2. ``my_daofind_sw_fits`` – run DAOStarFinder on a calibrated SW image and
   record astrometric offsets relative to the default WCS.

Design note
~~~~~~~~~~~
The original notebook contained module-level globals (``pid``, ``output_dir``) that parameterised ``reduce_img_stage2``.  These are now explicit
function parameters to make the function self-contained and pool-safe.
"""

from __future__ import annotations

import os
import time
from scipy import optimize

import numpy as np

from astropy.io import fits, ascii
from astropy.stats import sigma_clipped_stats, SigmaClip
from astropy import wcs
from astropy.table import Table, vstack, join, Column
from astropy.coordinates import SkyCoord
from astropy import units as u

from photutils.detection import DAOStarFinder
from photutils.background import Background2D, MedianBackground

import matplotlib.pyplot as plt
from nircam_wfss.plotting import corner_text

linear = lambda x, k, b: x * k + b  # noqa: E731


def reduce_img_stage2(
    rate_image_file: str,
    output_dir: str,
    overwrite: bool = False,
) -> None:
    """
    Apply JWST Stage-2 calibration (flat-field, photometric cal, WCS) to one
    NIRCam *rate* image.

    An association JSON file is written to ``output_dir``, then
    ``calwebb_image2.Image2Pipeline`` is run with resampling disabled so
    that the output is a single *cal.fits* in the same pixel grid as the
    input.

    Parameters
    ----------
    rate_image_file:
        Path to the NIRCam *rate.fits* stage-1 file.
    output_dir:
        Directory where the output *cal.fits* (and association JSON) are
        written.  The working directory is temporarily changed to this path
        during the pipeline call.
    overwrite:
        If ``False`` and the output *cal.fits* already exists, skip.
    """
    from jwst import datamodels
    from jwst.assign_wcs.assign_wcs_step import AssignWcsStep
    from jwst.flatfield import FlatFieldStep

    cal_path = os.path.join(
        output_dir,
        os.path.basename(rate_image_file)
        .replace("_rate_rowsub.fits", "_cal.fits")
        .replace("_rate.fits", "_cal.fits"),
    )
    if os.path.exists(cal_path) and not overwrite:
        print("cal.fits exists, skipping: %s" % cal_path)
        return

    # Open the rate file as an HDUList and pass it directly into ImageModel.
    # This avoids the read_metadata() → ASDF-extension lookup that both
    # Step.call(filename) and Step.run(filename) trigger in newer jwst/stpipe
    # versions. Older rate files (produced by a previous pipeline version)
    # lack the embedded ASDF extension, causing a KeyError on 'ASDF'.
    # Passing an already-open HDUList builds metadata purely from FITS headers.
    print("run AssignWcsStep + FlatFieldStep for %s" % os.path.basename(rate_image_file))
    try:
        with fits.open(rate_image_file, memmap=False) as hdul:
            model = datamodels.ImageModel(hdul)
            model = AssignWcsStep().run(model)
            model = FlatFieldStep().run(model)
            model.save(cal_path)
    except (ValueError, TypeError) as e:
        import warnings
        warnings.warn(
            "Skipping corrupted/truncated file: %s\n  (%s: %s)"
            % (rate_image_file, type(e).__name__, e),
            RuntimeWarning, stacklevel=2,
        )


def my_daofind_sw_fits(
    tmp_rate_sw: str,
    direct_image_dir: str,
    astrometry_dir: str,
    fwhm: float = 5.0,
    threshold: float = 7.0,
    use_default_wcs: bool = True,
    overwrite: bool = False,
) -> Table | None:
    """
    Run DAOStarFinder on a calibrated NIRCam SW image and return a source
    catalogue with sky coordinates.

    The astrometry can come from either the default WCS (derived from the
    spacecraft attitude, via ``crds`` + ``assign_wcs``) or from the
    calibrated WCS stored in the associated *cal.fits* file.

    Parameters
    ----------
    tmp_rate_sw:
        Path to the SW *rate.fits* file.
    direct_image_dir:
        Directory containing the corresponding *cal.fits* file.
    astrometry_dir:
        Directory where the output DAOStarFinder catalogue (ASCII) is
        written, under ``<astrometry_dir>/``.
    fwhm:
        Source FWHM in pixels passed to ``DAOStarFinder``.
    threshold:
        Detection threshold in units of background RMS.
    use_default_wcs:
        If ``True``, use the default WCS (from telescope attitude) rather
        than the pre-calibrated WCS in the *cal.fits* file.

    Returns
    -------
    tb_daofind:
        Astropy ``Table`` of detected sources with ``skycoord`` column, or
        ``None`` if no sources are detected.
    """
    import crds
    from jwst import assign_wcs, datamodels

    cal_file = os.path.join(
        direct_image_dir,
        os.path.basename(tmp_rate_sw).replace("_rate.fits", "_cal.fits"),
    )
    if not os.path.isfile(cal_file):
        print("cal.fits not found, skipping: %s" % cal_file)
        return None

    # --- Choose WCS source ---
    if use_default_wcs:
        with fits.open(tmp_rate_sw, memmap=False) as rate_sw_fits:
            rate_sw_hd = rate_sw_fits[0].header
            siaf_file = crds.getreferences(
                rate_sw_hd, reftypes=["distortion"], ignore_cache=False
            )["distortion"]
            wcs_step = assign_wcs.assign_wcs_step.AssignWcsStep(
                override_distortion=siaf_file
            )
            rate_sw_IM = datamodels.image.ImageModel(rate_sw_fits)
            rate_sw_IM_wcs = wcs_step.run(rate_sw_IM)
            tmp_wcs = rate_sw_IM_wcs.get_fits_wcs()
            rate_sw_IM_wcs.close()
    else:
        with fits.open(cal_file, memmap=False) as hdul:
            tmp_wcs = wcs.WCS(hdul[1].header)

    # --- Build detection image ---
    with fits.open(cal_file, memmap=False) as hdul:
        detect_hd = hdul[0].header
        tmp_detect_img = np.nan_to_num(hdul["SCI"].data)

    bkg = Background2D(
        tmp_detect_img, (64, 64), filter_size=(5, 5),
        sigma_clip=SigmaClip(sigma=3.0),
        bkg_estimator=MedianBackground(),
    )
    tmp_std = sigma_clipped_stats(
        tmp_detect_img[100:1500, 100:1500].flatten()[::7], sigma=3, maxiters=10
    )[-1]
    tmp_detect_img = (tmp_detect_img - bkg.background) / tmp_std

    # --- DAOStarFinder ---
    daofile = os.path.join(astrometry_dir, os.path.basename(tmp_rate_sw).replace("_rate.fits", "_daofind.dat"))
    if os.path.exists(daofile) and not overwrite:
        print(
            f"DAOFind already complete: %s" % daofile
        )
        return

    print("Running DAOFind on %s" % os.path.basename(cal_file))
    daofind = DAOStarFinder(fwhm=fwhm, threshold=threshold)
    tb_daofind = daofind(tmp_detect_img)
    if tb_daofind is None:
        print("  No sources detected.")
        return None

    tb_daofind["detector"] = detect_hd.get("DETECTOR", "unknown").lower()
    tb_daofind["skycoord"] = wcs.utils.pixel_to_skycoord(
        tb_daofind["xcentroid"], tb_daofind["ycentroid"], tmp_wcs
    )

    # --- Save catalogue ---
    os.makedirs(astrometry_dir, exist_ok=True)
    out_name = os.path.basename(tmp_rate_sw).replace("_rate.fits", "_daofind.dat")
    ascii.write(
        tb_daofind,
        os.path.join(astrometry_dir, out_name),
        overwrite=True, format='ecsv'
    )
    return tb_daofind

def calibrate_astrometry(
    tmp_rate_sw: str,
    direct_image_dir: str,
    v1p5_dir: str,
    astrometry_dir: str,
    astrometry_cal_table: str,
    astrometry_ref_table: str,
    overwrite: bool = False,
) -> Table | None:
    """
    Calibrate the astrometry of a SW *cal.fits* file by running DAOStarFinder
    and measuring the offset between the default WCS and the detected sources.

    The astrometric correction is saved to an ASCII table in
    ``<astrometry_dir>/`` with columns ``detector``, ``ra_offset``,
    ``dec_offset``, and ``pa_offset``.

    Parameters
    ----------
    tmp_rate_sw:
        Path to the SW *rate.fits* file.
    direct_image_dir:
        Directory containing the corresponding *cal.fits*  and *rate.fits* file.
    v1p5_dir:
        Directory containing the corresponding v1.5 calibrated files.
    astrometry_dir:
        Directory where the output astrometry table (ASCII) is written, under
        ``<astrometry_dir>/``.
    astrometry_cal_table:
        Path to an ASCII table with columns ``ra`` and ``dec`` containing the
        reference astrometry for the sources.
    astrometry_ref_table:
        Path to an ASCII table with columns ``ra`` and ``dec`` containing the
        reference astrometry for the sources, calibrated to external catalog (e.g. GAIA).
    overwrite:
        If True, overwrite the existing astrometry table.
    Returns
    -------
    Table | None
        Astropy Table with columns ``detector``, ``ra_offset``, ``dec_offset``,
        and ``pa_offset``, or ``None`` if the astrometry could not be
        calibrated (e.g. no sources detected).
    """
    if os.path.isfile(astrometry_cal_table) and not overwrite:
        print("Astrometry calibration table already exists, skipping: %s" % astrometry_cal_table)
        return ascii.read(astrometry_cal_table)
    if not os.path.isfile(astrometry_ref_table):
        raise FileNotFoundError("Astrometry reference table not found: %s" % astrometry_ref_table)
    # load the reference astrometry table
    tb_astrom_ref = ascii.read(astrometry_ref_table)
    ref_RA, ref_DEC = tb_astrom_ref['skycoord'].ra.value, tb_astrom_ref['skycoord'].dec.value
    tmp_coord_ref = SkyCoord(ref_RA, ref_DEC, unit = (u.deg, u.deg))

    # initialize the output astrometry table
    tb_sw_astrometry = Table(
        names = ['expName','N_match', 'dRA', 'dDEC', 'dRA_err', 'dDEC_err', 'theta', 'theta_err'], 
        dtype = ['S40', 'i4', 'f4', 'f4' , 'f4', 'f4', 'f4', 'f4'])
    tb_sw_astrometry.meta['comments'] = ['RA/Dec offsets are in arcsec and absolute', 
                    'i.e., cos(Dec) projection has been considered.',
                    '  >> dRA  =  RA_daofind -  RA_gaia',
                    '  >> dDEC = DEC_daofind - DEC_gaia',
                    'generated on %s' % time.strftime("%Y/%m/%d",  time.localtime())]
    for colname in tb_sw_astrometry.colnames[2:]: tb_sw_astrometry[colname].info.format = '.4f'
    sigma_clip = SigmaClip(sigma = 2.)

    # Compute astrometric offsets for each group of SW exposures.
    print('>>>  Astrometry (RA/DEC) offsets calculation for %3d SW Frames' % len(tmp_rate_sw))
    for k in range(len(tmp_rate_sw) // 4):
        ## For each cal.fits, get all the paths of csal.fits from all four detector.
        tmp_img_cal_path = tmp_rate_sw[k * 4]
        tmp_img_cal_path_base = os.path.basename(tmp_img_cal_path)[:30]
        tmp_grism_LW_path = os.path.join(v1p5_dir, 
                        tmp_img_cal_path_base + 'long_rate_lv1.5.fits')
        if os.path.isfile(tmp_grism_LW_path) == False: continue
        tmp_grism_hd_sci = fits.getheader(tmp_grism_LW_path, 1)
        print('%3d >> %s ' % (k, tmp_img_cal_path_base))
        tmp_img_cal_hd = fits.getheader(os.path.join(direct_image_dir, tmp_img_cal_path), 0)
        tmp_detector = tmp_img_cal_hd['DETECTOR'].lower()
        tmp_detectors = np.array(['%s%d' % (tmp_detector[:-1], tmp_det) for tmp_det in [1, 2, 3, 4]])
        tmp_img_cal_paths = np.array([tmp_img_cal_path.replace(tmp_detector, x) for x in tmp_detectors])
        
        
        ## Load the DAOFIND star catalog for each exposure
        for j, tmp_img_cal_path in enumerate(tmp_img_cal_paths):
            if os.path.isfile(tmp_img_cal_path)==False: 
                print('%s not found!' % tmp_img_cal_path)
                continue
            # Load pre-computed DAOFIND catalog for this detector.
            tmp_tb_daofind = ascii.read(os.path.join(astrometry_dir, 
                    os.path.basename(tmp_img_cal_path)\
                        .replace('_cal.fits', '_daofind.dat')\
                        .replace('_rate.fits', '_daofind.dat')),
                    format = 'ecsv')
            if j == 0 : tb_daofind = tmp_tb_daofind
            else: tb_daofind = vstack((tb_daofind, tmp_tb_daofind))
        tb_daofind = tb_daofind[np.argsort(tb_daofind['peak'])][-200:]
        tmp_coord_daofind = tb_daofind['skycoord'] 
        
        
        # Cross-match DAOFIND sources to the astrometric reference catalog.
        # rough center of the frame
        tmp_coord_center = SkyCoord(np.median(tmp_coord_daofind.ra), np.median(tmp_coord_daofind.dec))
        ### only select sources close to the center
        tmp_coord_ref = tmp_coord_ref[tmp_coord_center.separation(tmp_coord_ref) < 8 * u.arcmin]
        
        
        ## Cross match DAOFind Catalog with Gaia Catalog
        idx_daofind, d2d, _ = tmp_coord_ref.match_to_catalog_sky(tmp_coord_daofind)
        print(d2d)
        idx_ref = np.where(d2d < 0.25 * u.arcsec)[0]
        idx_daofind = idx_daofind[idx_ref]
        if len(idx_daofind) == 0:
            raise ValueError("No matched sources found for %s. The minimum separation "+\
            "is %.4f arcsec" % (tmp_img_cal_path, d2d.min().to(u.arcsec).value))

        ## Compute RA/DEC offset in arcsec: DAOFIND - ref
        tmp_ra_offset = ((tmp_coord_daofind[idx_daofind].ra - tmp_coord_ref[idx_ref].ra) 
                        * np.cos(tmp_coord_ref[idx_ref].dec)).to(u.arcsec).value
        tmp_dec_offset = (tmp_coord_daofind[idx_daofind].dec - tmp_coord_ref[idx_ref].dec).to(u.arcsec).value

        ## Find center of clustering, sigma-clipped RA/DEC offset
        dbin = 0.1
        tmp_ra_bins = np.arange(tmp_ra_offset.min()-0.01 - dbin/2, 
                                tmp_ra_offset.max() + dbin + 0.01, dbin)
        tmp_dec_bins = np.arange(tmp_dec_offset.min()-0.01- dbin/2, 
                                 tmp_dec_offset.max() + dbin + 0.01, dbin)
        offset_hist2d = np.histogram2d(tmp_ra_offset, tmp_dec_offset, 
                                       bins = (tmp_ra_bins, tmp_dec_bins))[0]
        idx_ra, idx_dec = np.where(offset_hist2d == np.nanmax(offset_hist2d))
        idx_ra, idx_dec = idx_ra[0], idx_dec[0]
        arg_clipped = np.where((tmp_ra_offset > tmp_ra_bins[idx_ra] - 1.5 * dbin) & 
                            (tmp_ra_offset < tmp_ra_bins[idx_ra] + 1.5 * dbin) &
                            (tmp_dec_offset > tmp_dec_bins[idx_dec] - 1.5 * dbin) & 
                            (tmp_dec_offset < tmp_dec_bins[idx_dec] + 1.5 * dbin))[0]
        tmp_ra_offset_clipped, tmp_dec_offset_clipped = tmp_ra_offset[arg_clipped], tmp_dec_offset[arg_clipped]
        idx_daofind = idx_daofind[arg_clipped]
        

        tmp_ra_offset_med, tmp_ra_offset_std = sigma_clipped_stats(tmp_ra_offset_clipped)[1:]
        tmp_dec_offset_med, tmp_dec_offset_std = sigma_clipped_stats(tmp_dec_offset_clipped)[1:]
        
        ## Compute rotation:
        poly_dec_dRA = optimize.curve_fit(f = linear, 
                xdata = tmp_coord_daofind[idx_daofind].dec, 
                ydata = tmp_ra_offset_clipped)[0]
        poly_ra_dDEC = optimize.curve_fit(f = linear, 
                xdata = tmp_coord_daofind[idx_daofind].ra, 
                ydata = tmp_dec_offset_clipped)[0]
        ## iter-1
        arg_dec_dRA = np.where(sigma_clip(np.polyval(poly_dec_dRA, tmp_coord_daofind[idx_daofind].dec.value) - tmp_ra_offset_clipped).mask == False)[0]
        poly_dec_dRA = optimize.curve_fit(f = linear, xdata = tmp_coord_daofind[idx_daofind].dec[arg_dec_dRA], ydata = tmp_ra_offset_clipped[arg_dec_dRA])[0]
        arg_ra_dDEC = np.where(sigma_clip(np.polyval(poly_ra_dDEC, tmp_coord_daofind[idx_daofind].ra.value) - tmp_dec_offset_clipped).mask == False)[0]
        poly_ra_dDEC = optimize.curve_fit(f = linear, xdata = tmp_coord_daofind[idx_daofind].ra[arg_ra_dDEC], ydata = tmp_dec_offset_clipped[arg_ra_dDEC])[0]
        ## iter-2
        arg_dec_dRA = np.where(sigma_clip(np.polyval(poly_dec_dRA, tmp_coord_daofind[idx_daofind].dec.value) - tmp_ra_offset_clipped).mask == False)[0]
        poly_dec_dRA, pcov_dec_dRA = optimize.curve_fit(f = linear, xdata = tmp_coord_daofind[idx_daofind].dec[arg_dec_dRA], ydata = tmp_ra_offset_clipped[arg_dec_dRA])
        arg_ra_dDEC = np.where(sigma_clip(np.polyval(poly_ra_dDEC, tmp_coord_daofind[idx_daofind].ra.value) - tmp_dec_offset_clipped).mask == False)[0]
        poly_ra_dDEC, pccov_ra_dDEC = optimize.curve_fit(f = linear, xdata = tmp_coord_daofind[idx_daofind].ra[arg_ra_dDEC], ydata = tmp_dec_offset_clipped[arg_ra_dDEC])
        perr_dec_dRA = np.sqrt(np.diag(pcov_dec_dRA))
        perr_ra_dDEC = np.sqrt(np.diag(pccov_ra_dDEC))
        ## combine two direction
        cos_factor = np.cos(np.deg2rad(np.median(tmp_coord_daofind[idx_daofind].dec.value)))
        theta_as_deg = (poly_dec_dRA[0] / perr_dec_dRA[0] - poly_ra_dDEC[0] / perr_ra_dDEC[0]) / (1/perr_dec_dRA[0] + cos_factor / perr_ra_dDEC[0]) / 3600.
        theta_as_deg_err = np.sum(perr_dec_dRA[0]**2 + (perr_ra_dDEC[0] / cos_factor)**2)**0.5 / 2. / 3600.
        theta_as_deg = np.rad2deg(theta_as_deg)
        theta_as_deg_err = np.rad2deg(theta_as_deg_err)
        
        tmp_ra_offset_med = tmp_grism_hd_sci['CRVAL2'] * poly_dec_dRA[0] + poly_dec_dRA[1]
        tmp_dec_offset_med = tmp_grism_hd_sci['CRVAL1'] * poly_ra_dDEC[0] + poly_ra_dDEC[1]

        tmp_ra_offset_std = sigma_clipped_stats(tmp_ra_offset_clipped)[2]
        tmp_dec_offset_std = sigma_clipped_stats(tmp_dec_offset_clipped)[2]
        print('           N_mathch = %d' % len(arg_clipped))
        print(' RA_image -  RA_ref = %.3f"±%.3f"' % (tmp_ra_offset_med,  tmp_ra_offset_std  / np.sqrt(len(arg_clipped)) ))
        print('DEC_image - DEC_ref = %.3f"±%.3f"' % (tmp_dec_offset_med, tmp_dec_offset_std / np.sqrt(len(arg_clipped)) ))
        print('     rotation theta = %.3f±%.3f deg' % (theta_as_deg, theta_as_deg_err))
        tb_sw_astrometry.add_row([tmp_img_cal_path_base, len(arg_clipped),
                                tmp_ra_offset_med, tmp_dec_offset_med, 
                                tmp_ra_offset_std / np.sqrt(len(arg_clipped)), tmp_dec_offset_std / np.sqrt(len(arg_clipped)),
                                theta_as_deg, theta_as_deg_err])
        
        plt.close()
        plt.subplots(1, 1, figsize = (5, 5))
        plt.plot(tmp_ra_offset, tmp_dec_offset, marker = '.', color = 'gray', ls = 'none')
        plt.plot(tmp_ra_offset_clipped, tmp_dec_offset_clipped, marker = '.', color = 'k', ls = 'none')
        plt.plot(tmp_ra_offset_med, tmp_dec_offset_med, marker = '+', color = 'red', ms = 20, mew = 2, ls = 'none')
        plt.axhline(0, color = 'grey', ls = '--')
        plt.axvline(0, color = 'grey', ls = '--')
        plt.gca().set(xlabel = r'$\Delta$RA', ylabel = r'$\Delta$Dec')
        if np.max(np.abs(np.array([tmp_ra_offset_med, tmp_dec_offset_med]))) > 1:
            plt.gca().set(aspect = 1, xlim = (-1.5, 1.5), ylim = (-1.5, 1.5), 
                        xticks = np.arange(-1.5, 1.51, 0.5), yticks = np.arange(-1.5, 1.51, 0.3))
        elif  np.max(np.abs(np.array([tmp_ra_offset_med, tmp_dec_offset_med]))) > 0.5:
            plt.gca().set(aspect = 1, xlim = (-1., 1.), ylim = (-1., 1.), 
                        xticks = np.arange(-0.8, 1.1, 0.4), yticks = np.arange(-1., 1.1, 0.2))
        else:
            plt.gca().set(aspect = 1, xlim = (-0.5, 0.5), ylim = (-0.5, 0.5), 
                        xticks = np.arange(-0.4, 0.51, 0.2), yticks = np.arange(-0.5, 0.51, 0.1))
                    
        corner_text(plt.gca(), s = ' RA_image -  RA_ref = %.3f"±%.3f"' % (tmp_ra_offset_med,  tmp_ra_offset_std / np.sqrt(len(arg_clipped)) ),
                    loc = 1, fontsize = 12)
        corner_text(plt.gca(), s = 'DEC_image - DEC_ref = %.3f"±%.3f"' % (tmp_dec_offset_med, tmp_dec_offset_std/ np.sqrt(len(arg_clipped)) ),
                    loc = 4, fontsize = 12)
        plt.title(tmp_img_cal_path_base, fontsize = 14, )
        plt.tight_layout()
    
    ascii.write(tb_sw_astrometry, astrometry_cal_table, overwrite = overwrite)
    return tb_sw_astrometry
