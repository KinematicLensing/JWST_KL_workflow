# Measuring Kinematic Lensing with JWST FRESCO+JADES Data

This repository archives the workflow of a preliminary exploration of kinematic lensing (KL) measurement using JWST NIRCam grism data, FRESCO, and JWST imaging data (JADES). 
To run the KL analysis, you need to install the kl-tools and download the FRESCO+JADES rate data. 

## 1. Python Environment
This KL analysis is build on top of two kinds of python environment: `jwst` for JWST data reduction, and `kltools` for KL fitting. 

### 1.1 `kl-tools` Installation
The environment of `kltools` is recorded in `environments/kltools_environment.yml`. After creating this Python env, please install [kltools](https://github.com/sweverett/kl-tools/tree/tng_grism) (NOTE: use the `tng_grism` branch) by

``` 
conda activate ${YOUR_KLTOOLS_ENV}
git clone https://github.com/sweverett/kl-tools/tree/tng_grism 
cd kl_tools/grism_modules/src
python setup.py install
cd ../../../
python setup.py install
``` 

### 1.2 `jwst` Environment
Create the `jwst` data reduction environment from `environments/jwst_environment.yml`

## 2. JWST Data Download

### 2.1 Download FRESCO & CONGRESS Data

See the notebook [Download_FRESCO_CONGRESS_Data.ipynb](./notebook/Download_FRESCO_CONGRESS_Data.ipynb) for FRESCO and CONGRESS grism and direct imaging data downloads. The datasets are saved to:
- FRESCO grism: `/xdisk/timeifler/jiachuanxu/jwst/fresco/grism_F444W` 
- FRESCO direct image: `/xdisk/timeifler/jiachuanxu/jwst/fresco/direct_imaging`
- CONGRESS grism: `/xdisk/timeifler/jiachuanxu/jwst/congress/grism_F356W`
- CONGRESS direct image: `/xdisk/timeifler/jiachuanxu/jwst/congress/direct_imaging`


### 2.2 Download JADES Data

The JADES DR5 imaging mosaics can be downloaded from [here](https://slate.ucsc.edu/~brant/jades-dr5/). We will be mostly interested in F444W band since that's the same band as the grism image, and the morphology of galaxies in NIR bands are more smooth. But we download F090W/F200W/F444W mosaics to build a pseudo-RGB stamp image for the galaxies as sanity check. 

The JADES multiband drizzled mosaics are downloaded to `/xdisk/timeifler/jiachuanxu/jwst/jades/mosaics/hlsp_jades_jwst_nircam_${FIELD}-deep_${BAND}_v5.0_drz.fits`. 

<!--
Since the JADES mosaic files are large, 2D multi-band cutouts can also be requested per object given the RA and DEC of the galaxy. See [Request_JADES_Cutouts.ipynb] for an example.~ Not working, need debug
-->

### 2.3 Download JADES Catalog

The new DR5 catalog can be downloaded [here](https://slate.ucsc.edu/~brant/jades-dr5/), which includes better detection, deblending, and phtometry. Consider using this catalog to decide if an object is blended, and select isoloated objects only. (The `PARENT_ID` in the `FLAG` extension, and also search for nearest objects, see if their Kron radius overlap)

**TODO** After some trial and error, the JADES DR5 deblending seems not quite reliable. Consider trying some other deblending methods, like SCARLET, in future, but we'll stick with JADES DR5 catalog for now. 

For catalog documentation, see the appendix of [Robertson et al. (2026)](https://arxiv.org/pdf/2601.15956)


## 3. JWST Grism Data Reduction

We will run our own grism calibration and mosaic, which is modified from [Fengwu Sun's JWST NIRCam WFSS data reduction pipeline](https://github.com/fengwusun/nircam_grism). The modified version is [xxx](xxx). Major changes are:
- Skip sensitivity calibration during data reduction and keep the `DN/s` unit, because we forward model the sensitivity.
- We coadd the grism image grouped by roll angles / dispersion directions. 
- blabla

To reduce the grism data, we first get a catalog of line emitting galaxies, do a preliminary target selection to select galaxies that are useful for KL measurement. Then, we will calculate the grism dispersion trace of each object, and extract their 2D grism spectra. 2D spectra are coadded grouped by dispersion angle, and a continuum is subtracted by mediam filter. We will build 2D cutouts including their brandband image and emission line-only 2D grism. 

### 3.1 Catalog Adjustment

The low-z Paschen and Brackett emitters in the original Fengwu's catalog may have artifacts in centroid determination and may contain duplicated objects. Therefore, we need to match Fengwu's catalog to JADES public catalog (DR5) and adjust the RA/DEC and OBJID.

The original Fengwu's catalogues are 
- `v094_gds_fresco_line_list_low-ground-z_Pa_Br.fits`, 
- `v091_gdn_congress_line_list_v2_low-ground-z_Pa_Br.fits`, 
- `v091_gdn_fresco_line_list_low-ground-z_Pa_Br.fits`.

After removing duplicated objects and unmatched objects, the number of galaxies are:
- FRESCO GOODS-S: 296 -> 281
- FRESCO GOODS-N: 260 -> 251
- CONGRESS GOODS-N:  345 -> 336

The adjusted catalogs are saved to `/xdisk/timeifler/jiachuanxu/jwst/fengwu_catalog`:  
- `v1_gds_fresco_line_list_low-ground-z_Pa_Br_jades_dr5.fits`
- `v1_gdn_fresco_line_list_low-ground-z_Pa_Br_jades_dr5.fits`
- `v1_gdn_congress_line_list_low-ground-z_Pa_Br_jades_dr5.fits`. 


### 3.2 Imaging Target Selection

We don't need to reduce the grism spectra for all galaxies because, e.g., galaxies in merging or blending systems are not useful. We do the first round of **target selection** based on the broadband image of those galaxies to save computational resources --- but feel free to extract the 2D spectrum of discarded objects for your interest. 

KL sample first needs to pass standard weak lensing source sample selectrion criteria:

- No blending/merging galaxies: for now, we require that the Kron ellipse of the 5th nearest objects do not intersect the Kron ellipse of the primary target. Otherwise, if there's overlap, we require that the total flux from contaminant objects is no more than 10 percent of the primary object flux.  
- No bleeding (visual examination)
- Star-galaxy separation shows it's an extended galaxy, not point source like star or AGN
- Not on the edge of detector (visual examination)
- spatial resolution factor R > 0.4 (star-galaxy separation and spatial resolution factor is roughly the same thing...)

We down-select the v1 catalog produced in [Sect. 3.1](#31-catalog-adjustment) based on these criteria (see notebook [`Preliminary_Target_Selection.ipynb`](./notebook/Preliminary_Target_Selection.ipynb)), and save the catalog to `v2` version. 

A rough number of decrease as we apply some of the criteria
- Spatial resolution factor R > 0.4:
    - FRESCO GOODS-N: 251 -> 251
    - FRESCO GOODS-S: 281 -> 281
    - CONGRESS GOODS-N: 336 -> 336
- No blending:
    - FRESCO GOODS-N: 251 -> 218
    - FRESCO GOODS-S: 281 -> 219
    - CONGRESS GOODS-N: 336 -> 280

So all galaxies pass the WL size criteria, and roughly 80 percent of galaxies survived after blending cut. Note that we will re-iterate on the blending detection in future when we have better blending detection algo. 


### 3.3 Grism Reduction

For the objects in v2 catalog, we extract their 2D grism spectrum by running the customized NIRCam/WFSS pipeline [xxx?](xxx?)


## 4. Target Selection



## Food for thought

- Galaxy morphology breakdown: Since we are in an era of JWST + AI, why not classify galaxies based on their multi-band difraction-limited image? There must be some tool to return binary or discrete flag about a galaxy includes strong bar/spiral arms ([ZooBot](https://github.com/mwalmsley/zoobot)). Some code should be able to mask out or remove galactic bulge, bar, and spiral arm ([SpArcFiRe](https://arxiv.org/pdf/1707.02021)).
- The GOODS fields are too small, such that survey geometry might produce strong enough leakage into B-mode. Think about how to mitigate this
- PSF residual: we can't assume the webbpsf is perfect. Another PSF calibration method is used in [a JWST WL paper](https://arxiv.org/pdf/2304.02054)



