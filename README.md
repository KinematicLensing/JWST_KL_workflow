# Measuring Kinematic Lensing with JWST FRESCO+JADES Data

This repository archives the workflow of a priliminary exploration of kinematic lensing (KL) measurement using JWST NIRCam grism data, FRESCO, and JWST imaging data (JADES). 
To run the KL analysis, you need to install the kl-tools and download the FRESCO+JADES rate data. 


## 1. `kl-tools` Installation

## 2. JWST Data Download

### 2.1 Download FRESCO Data

See the notebook [Download_FRESCO_Data.ipynb] for FRESCO grism and direct imaging data download. Note that the code 

The FRESCO grism data is saved to `/xdisk/timeifler/jiachuanxu/jwst/fresco/grism_F444W` and the direct imaging data are saved to `/xdisk/timeifler/jiachuanxu/jwst/fresco/direct_imaging`

### 2.2 Download CONGRESS Data


### 2.3 Download JADES Data

The JADES imaging mosaics can be downloaded from [MAST](https://archive.stsci.edu/hlsp/jades#section-268de08a-1ff5-430e-adfe-846e6b933f3b). We will be mostly interested in F444W band since that's the same band as the grism image, and the morphology of galaxies in NIR bands are more smooth. But we download F090W/F200W/F444W mosaics to build a psedup-RGB stamp image for the galaxies as sanity check. 

The JADES multiband mosaics are downloaded to `/xdisk/timeifler/jiachuanxu/jwst/jades/mosaics/hlsp_jades_jwst_nircam_${FIELD}-deep_${BAND}_v2.0_drz.fits`. 

~Since the JADES mosaic files are large, 2D multi-band cutouts can also be requested per object given the RA and DEC of the galaxy. See [Request_JADES_Cutouts.ipynb] for an example.~ Not working, need debug

Note that the DR2 contains mosaics of GOODS-S **including** the F182M, F210M, and F444W from FRESCO, and the DR3 contains mosaics of GOODS-N **including** the F182M, F210M, and F444W from FRESCO. The original JADES proposal does not cover the full FoV of FRESCO. 

### 2.4 Download JADES Catalog

The JADES contains a series of [catalogs](https://archive.stsci.edu/hlsp/jades#section-268de08a-1ff5-430e-adfe-846e6b933f3b) of photometry and size measurements, see the [documentation for catalog extensions](https://archive.stsci.edu/hlsps/jades/hlsp_jades_jwst_nircam_goods-s-deep_photometry_v2.0_catalog-ext-readme.pdf).

NOTE: there's new DR5 catalog [here](https://slate.ucsc.edu/~brant/jades-dr5/) which includes better detection, deblending, and phtometry. Consider using this catalog to decide if an object is blended, and select isoloated objects only. (The `PARENT_ID` in the `FLAG` extension, and also search for nearest objects, see if their Kron radius overlap)

**TODO** After some trial and error, the JADES DR5 deblending seems not quite reliable. Consider trying some other deblending methods, like SCARLET, in future, but we'll stick with JADES DR5 catalog for now. 

For catalog documentation, see the appendix of [Robertson et al. (2026)](https://arxiv.org/pdf/2601.15956)


## 3. JWST Grism Data Reduction

We will run our own grism calibration and mosaic, which is modified from [Fengwu Sun's JWST NIRCam WFSS data reduction pipeline](https://github.com/fengwusun/nircam_grism). The modified version is [xxxxx.ipynb]. Major changes are:
- blabla
- blabla
- blabla

### 3.1 Catalog Adjustment

The low-z Paschen and Brackett emitters in the original Fengwu's catalog may have artifacts in centroid determination and may contain duplicated objects. Therefore, we need to match Fengwu's catalog to JADES public catalog (DR5) and adjust the RA/DEC and OBJID.

The original Fengwu's catalogues are `v094_gds_fresco_line_list_low-ground-z_Pa_Br.fits`, `v091_gdn_congress_line_list_v2_low-ground-z_Pa_Br.fits`, `v091_gdn_fresco_line_list_low-ground-z_Pa_Br.fits`.

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

We down-select the v1 catalog produced in 3.1 based on these criteria (see notebook `Image_Target_Selection.ipynb`), and save the catalog to `v2` version. 

A rough number of decrease as we apply some of the criteria
- Spatial resolution factor R > 0.4:
    - FRESCO GOODS-N: 
    - FRESCO GOODS-S:
    - CONGRESS GOODS-N:
- No blending:
    - FRESCO GOODS-N: 
    - FRESCO GOODS-S:
    - CONGRESS GOODS-N:


### 3.3 Grism Reduction

## 4. Target Selection






