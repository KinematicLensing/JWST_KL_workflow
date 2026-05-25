# NIRCam WFSS Grism Extraction Pipeline

A Python package for extracting 2D and 1D spectra from JWST NIRCam Wide-Field Slitless Spectroscopy (WFSS) grism observations.

**Acknowledgement.** This pipeline is derived from F. Sun's (CfA) original notebook-based implementation at [fengwusun/nircam_grism](https://github.com/fengwusun/nircam_grism). For a pedagogical introduction to the instrument mode and the extraction methodology, refer to that repository first. The code here refactors the original notebook into a production-ready Python package with a configuration-driven CLI, parallel processing, and clean module boundaries.

---

## Table of Contents

1. [NIRCam WFSS Background](#1-nircam-wfss-background)
2. [Installation](#2-installation)
3. [Quick Start](#3-quick-start)
4. [Pipeline Workflow](#4-pipeline-workflow)
5. [Configuration Reference](#5-configuration-reference)
6. [Package Layout](#6-package-layout)

---

## 1. NIRCam WFSS Background

### Instrument modes

NIRCam WFSS places a grism in the pupil wheel of the Long-Wavelength (LW) channel. Two grisms are available in each of the two modules (A and B), and both modules observe simultaneously, giving up to four independent grism datasets per visit:

| Label | Module | Dispersion axis | Wavelength range |
|-------|--------|-----------------|------------------|
| ModA-GRISMR | A | Row (horizontal) | 2.4–5.0 µm (F444W) |
| ModA-GRISMC | A | Column (vertical) | 2.4–5.0 µm (F444W) |
| ModB-GRISMR | B | Row (horizontal) | 2.4–5.0 µm (F444W) |
| ModB-GRISMC | B | Column (vertical) | 2.4–5.0 µm (F444W) |

The usable wavelength range shifts with the blocking filter: F356W covers roughly 3.1–4.0 µm, while F444W covers 3.8–5.2 µm. The LW pixel scale is 0.063 arcsec/pixel.

Accompanying NIRCam SW direct-imaging exposures (taken simultaneously via the dichroic) are used solely for astrometric calibration; they are not co-added into the spectral products.

### Pick-off mirror (POM) vignetting

The GRISM pick-off mirror does not fully overlap with the LW field of view. Sources near the edges of the detector, or whose spectra cross the POM boundary, receive only partial spectral coverage or are entirely obscured. The pipeline tracks POM transmission per source per frame using pre-computed transmission maps.

### How 2D spectral extraction works intuitively

Unlike a slit spectrograph, a grism disperses *every* source in the field simultaneously onto the same detector. The key observation is that the dispersion solution is a deterministic polynomial function of a source's detector position: given where a source sits on the detector (x₀, y₀), the pipeline can predict, for every wavelength λ, exactly which pixel (x₀ + dx(λ), y₀ + dy(λ)) that wavelength falls on.

The extraction therefore proceeds source by source from the same set of grism exposures:

1. **Locate the source.** Use the source's RA/Dec and the per-exposure WCS to project the source onto detector coordinates (x₀, y₀).
2. **Compute the spectral trace.** The dispersion polynomial maps wavelength → column offset dx; a trace polynomial maps (dx, position) → row offset dy. Together these define the curved path that the spectrum follows across the detector.
3. **Cut out the 2D stamp.** A rectangular aperture of half-width `aperture_pix` pixels is extracted perpendicular to the dispersion direction, centred on the trace, for each wavelength column.
4. **Assign wavelengths.** Each column of the 2D stamp is labelled with its corresponding wavelength via the inverted dispersion relation, producing a uniformly-sampled 2D spectrum in (wavelength, cross-dispersion) space.
5. **Co-add exposures.** Multiple exposures (from different visits or both grisms) are sigma-clipped and weighted-averaged to produce the final 2D spectrum.

For emission-line science the pipeline additionally separates continuum and emission-line flux (Stage 2e) before extraction, so the EMLINE extension of each grism exposure contains only line flux, suppressing continuum contamination from overlapping sources.

---

## 2. Installation

### Prerequisites

The pipeline depends on the JWST calibration pipeline (`jwst`), CRDS, and several standard astronomy packages. The recommended way to set up the environment is with conda:

```bash
conda create -n nircam_wfss python=3.11
conda activate nircam_wfss
pip install jwst crds astropy scipy matplotlib drizzle stpsf photutils
```

You also need to configure the CRDS environment before running:

```bash
export CRDS_SERVER_URL=https://jwst-crds.stsci.edu
export CRDS_PATH=/your/local/crds/cache      # path where CRDS files are cached
```

### Install this package

```bash
cd /path/to/JWST_KL/nircam_wfss
pip install -e .
```

This installs the `nircam_wfss` package and registers the `run-pipeline` CLI entry point.

### Calibration support data

The pipeline requires a set of pre-computed calibration files in `data/FSun_cal/` and `data/GRISM_NIRCAM/` (included in this repository):

| Directory | Contents |
|-----------|----------|
| `data/FSun_cal/` | Dispersion solution coefficients, sensitivity curves, median super-sky backgrounds, POM spectral coverage maps, astrometry reference catalogs |
| `data/GRISM_NIRCAM/` | POM transmission maps for modules A and B |
| `data/jwst_pipeline_config/` | JWST pipeline configuration files |

---

## 3. Quick Start

1. **Copy and edit a config file.** Start from an example in `configs/`:

   ```bash
   cp configs/PID1895_FRESCO_GDN_F444W.yaml configs/my_survey.yaml
   ```

   Edit the paths and parameters for your dataset (see [Section 5](#5-configuration-reference)).

2. **Run the pipeline:**

   ```bash
   run-pipeline configs/my_survey.yaml
   # or equivalently:
   python src/run_pipeline.py configs/my_survey.yaml
   ```

3. **Outputs** land in the directories you specified:
   - `calibrated_dir/lv1p5/` — WCS-assigned, flat-fielded, background-subtracted grism exposures
   - `calibrated_dir/plots/` — diagnostic PNGs for each processing step
   - `calibrated_dir/POM_catalog/` — per-frame source position catalogs with POM flags
   - `extract_dir/` — per-source 2D spectrum FITS files, 1D spectrum tables, and diagnostic PDFs

---

## 4. Pipeline Workflow

The pipeline runs eleven stages in sequence. Stages 2a–2e pre-process every grism exposure. Stages 3–5 prepare astrometry and source catalogs. Stages 6–8 extract spectra. Computationally expensive stages use a multiprocessing pool with `n_procs` workers.

```
rate.fits (JWST Stage-1 output)
        │
        ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│ Stage 2a  Assign WCS + flat-field         rate.fits → lv1.5.fits            │
│ Stage 2b  Build median super-sky          lv1.5.fits → bkg_*.fits           │
│ Stage 2c  Background subtraction          lv1.5.fits updated in place       │
│ Stage 2d  Hot-pixel rejection             lv1.5.fits updated in place       │
│ Stage 2e  Continuum subtraction           lv1.5.fits + EMLINE extension     │
└─────────────────────────────────────────────────────────────────────────────┘
        │
        ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│ Stage 3   SW direct image reduction       rate.fits → cal.fits (SW only)    │
│ Stage 4   Astrometry calibration          cal.fits → dRA/dDec table         │
│ Stage 5   Build POM-applied catalogs      source cat → per-frame catalogs   │
└─────────────────────────────────────────────────────────────────────────────┘
        │
        ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│ Stage 6   Extract 2D spectra              → spec_2d_*.fits                  │
│ Stage 7   Extract 1D spectra              → spec_1d_*.fits + diagnostic PDF │
│ Stage 8   Extract 2D emission-line stamps → emline_2d_*.fits (optional)     │
└─────────────────────────────────────────────────────────────────────────────┘
```

### Stage 2a — Assign WCS and flat-field

**Goal:** Convert each raw Stage-1 `rate.fits` grism file into a `lv1.5.fits` file that carries a valid WCS and has been flat-fielded.

**How:** The JWST `AssignWcsStep` and `FlatFieldStep` are applied. Because the JWST pipeline's WCS assignment rejects GRISM-mode exposures, each file is temporarily re-labelled as a direct-imaging exposure (`EXP_TYPE = NRC_IMAGE`) so that the pipeline fetches the correct distortion and flat-field reference files from CRDS, then restores the original metadata. The output `lv1.5.fits` contains SCI, ERR, and DQ extensions with WCS keywords in the header.

**Controlled by:** `overwrite`

---

### Stage 2b — Build median super-sky backgrounds

**Goal:** Create a median-stacked "super-sky" background image for each (filter, module, pupil) combination.

**How:** All `lv1.5.fits` files sharing the same filter/module/pupil are loaded. A 2σ-clipped median is computed pixel-by-pixel, suppressing astrophysical sources that would contaminate a simple mean. For GRISMC (column-dispersed), a special post-processing step handles the coronagraph transition region. The resulting backgrounds are saved as `median_bkg_<filter>_<module>_<pupil>.fits` in `cali_support_dir`. If a pre-computed background already exists (from `cali_support_dir`), this stage is skipped.

**Controlled by:** `cali_support_dir`

---

### Stage 2c — Background subtraction and 1/f noise removal

**Goal:** Subtract the sky background and remove detector 1/f correlated read noise from each exposure.

**How (per exposure):**

1. The median super-sky image is scaled to match each exposure's sky level and subtracted.
2. A 2D residual background is estimated using SExtractor-style segmentation, masking detected sources, and fitting a smooth background to the remaining pixels.
3. 1/f noise (horizontal or vertical striping) is removed by computing the per-row or per-column median of source-masked pixels and subtracting it. The pipeline supports both single-channel (full-row) and 4-amplifier-channel (512-pixel segment) subtraction modes.

The `lv1.5.fits` SCI array is updated in place. An optional 4-panel diagnostic PNG is written to `plot_dir`.

**Controlled by:** `cali_support_dir`

---

### Stage 2d — Hot-pixel rejection

**Goal:** Identify and flag isolated hot pixels that were not caught by the Stage-1 jump detection.

**How:** A roll-difference method is applied: each pixel is compared to its two neighbours along the dispersion direction. If the difference exceeds `sigma_hot` times the local noise, the pixel is flagged in the DQ extension and its SCI and ERR values are replaced with NaN. This approach is sensitive to isolated spikes rather than extended structures, so real astronomical emission lines are preserved.

**Controlled by:** `sigma_hot`

---

### Stage 2e — Continuum subtraction and emission-line separation

**Goal:** Separate the broadband stellar/galaxy continuum from narrow emission lines in each grism exposure, enabling cleaner emission-line-only extraction in Stage 6.

**How:** A two-pass median-filter approach is used along the dispersion direction:

1. **Pass 1:** A wide median filter with a central hole (to avoid self-subtraction) is applied along each row (GRISMR) or column (GRISMC). This isolates the smooth continuum, and subtracting it leaves the emission lines. Residual 1/f noise is also removed at this step.
2. **Pass 2:** High signal-to-noise emission pixels from Pass 1 are masked, and the median filter is re-run on the masked data for a cleaner continuum estimate.

Two new FITS extensions are appended to the `lv1.5.fits` file: `EMLINE` (line-only flux) and `CONT` (continuum model). Stage 6 uses the `EMLINE` extension when extracting emission-line cutouts in Stage 8.

---

### Stage 3 — SW direct image Stage-2 reduction

**Goal:** Apply JWST Stage-2 calibration to the SW direct-imaging exposures so they can be used for astrometric calibration.

**How:** Each SW `rate.fits` file is processed through `AssignWcsStep` and `FlatFieldStep` (no resampling) to produce a `cal.fits` file with a WCS and photometric calibration. This step is skipped if no SW rate files are found.

**Controlled by:** `direct_image_dir`, `direct_image_filename_fmt`

---

### Stage 4 — Astrometry calibration

**Goal:** Measure and record the astrometric offset (ΔRA, ΔDec, rotation θ) between the JWST pointing model and an external astrometric reference (e.g., Gaia via HST).

**How:** If a pre-computed astrometry table already exists at `astrometry_cal_table`, it is loaded directly (the stage is effectively skipped). Otherwise, DAOStarFinder is run on each SW `cal.fits` to detect point sources, and the resulting catalogs are cross-matched to the external reference catalog at `astrometry_ref_table`. Sigma-clipped statistics yield per-exposure ΔRA, ΔDec, and θ values, which are saved to `astrometry_cal_table` for reuse. These corrections are applied to every source position in Stages 5 and 6.

**Controlled by:** `astrometry_cal_table`, `astrometry_ref_table`, `direct_image_dir`

---

### Stage 5 — Build POM-applied catalogs

**Goal:** For each grism exposure, determine which sources from the input catalog are observable (not vignetting by the POM) and compute their exact detector positions with astrometric corrections applied.

**How:** For each source in the input catalog, its RA/Dec is projected into detector coordinates for every `lv1.5.fits` frame, applying the per-exposure astrometric corrections from Stage 4. The POM transmission maps are consulted to determine whether the source itself, and what fraction of its spectrum, falls within the unvignetted region. Sources below the POM coverage threshold are excluded. Per-frame catalogs are saved in `POM_catalog_dir` as FITS tables, and two dictionaries are returned in memory: one mapping each source ID to the list of frame indices where it is observable, and one mapping it to the corresponding catalog file paths.

**Controlled by:** `source_catalog_path`, `cali_support_dir`, `default_POM_trans_dir`

---

### Stage 6 — Extract 2D grism spectra

**Goal:** Produce a per-source 2D spectrum FITS file containing contributions from all observable grism exposures, co-added by module/pupil and then globally.

**How (per source, per exposure):**

1. The source position and astrometric corrections are read from the POM catalog for this frame.
2. The grism configuration object (`GrismConf`) is used to compute the spectral trace: the polynomial dispersion relation maps wavelength → dx pixel offset, and the trace polynomial gives the cross-dispersion offset dy. Together they define the curved path the spectrum follows across the detector.
3. A rectangular 2D stamp of width 2 × `aperture_pix` pixels is cut out perpendicular to the trace at each wavelength column, using sub-pixel interpolation to follow the curved trace.
4. For GRISMC (vertical dispersion), the stamp is transposed so that the wavelength axis always runs horizontally in the output.
5. The resulting stamps from all exposures are assembled into a single multi-extension FITS file per source. Extensions are labelled `SPEC2D-N`, `WHT2D-N`, `DQ2D-N`, `LINE2D-N` (emission-line channel), and `WAVE` (wavelength array). Co-additions are formed separately per module/pupil (to allow grism-cross-checks) and then across all combinations, using sigma-clipped inverse-variance weighting. Units are converted according to `bunit_spec2d`.

**Controlled by:** `aperture_pix`, `bunit_spec2d`, `overwrite_spec2d`, `source_catalog_path`

---

### Stage 7 — Extract 1D spectra

**Goal:** Collapse the co-added 2D spectrum into a 1D spectrum using optimal (Horne) or boxcar extraction, and produce a diagnostic PDF for quality assessment.

**How:** The all-module-pupil co-added 2D spectrum from Stage 6 is read. The spatial profile (cross-dispersion profile) is estimated either from a direct-image cutout (read from `image_mosaic_dir`) or from the collapsed 2D spectrum itself. A Gaussian is fitted to the profile to determine the source centroid and width. Two 1D extractions are performed: a boxcar (aperture of ±3σ) and an optimal (profile-weighted) extraction. A diagnostic PDF is saved to `plot_dir` showing the 2D spectrum, the spatial profile, and the extracted 1D spectrum alongside an RGB direct-image thumbnail.

**Controlled by:** `image_mosaic_dir`, `image_mosaic_filename_fmt`, `image_mosaic_rgb_bands`, `image_mosaic_field`

---

### Stage 8 — Extract 2D emission-line cutouts

**Goal:** Produce a small 2D cutout centred on a known emission line for each source, co-added across all grism exposures, for morphological or kinematic analysis.

**How:** This stage requires `z_grism` (spectroscopic redshift) and `name_line_exp` (emission-line name) columns in the source catalog. The predicted detector position of the emission line is computed using the dispersion solution and the source's redshift. A cutout of size `cutout_size_simple` or `cutout_size_drizzle` pixels is extracted from the `EMLINE` extension of each `lv1.5.fits` exposure. Two co-addition modes are available:

- **`simple`**: Cutouts are stacked directly in detector frame without re-sampling, preserving the native pixel scale. Each frame's position angle (from `GS_V3_PA` header) is recorded in the output header for downstream derotation.
- **`drizzle`**: Each cutout is re-sampled onto a common equatorial grid using the Drizzle algorithm, with configurable output pixel scale (`finalscale_drizzle`) and drop-size fraction (`pixfrac_drizzle`). This produces north-up, east-left co-added stamps at sub-pixel resolution.

**Controlled by:** `coadd_method`, `cutout_size_simple`, `cutout_size_drizzle`, `finalscale_drizzle`, `pixfrac_drizzle`, `psf_oversample`

---

## 5. Configuration Reference

The YAML configuration file controls every aspect of the pipeline. Unrecognised keys are silently ignored. All paths may be absolute or relative to the working directory.

### Programme metadata

| Key | Type | Description | Used in Stage |
|-----|------|-------------|---------------|
| `pid` | int | JWST Programme ID (cosmetic; used in association files) | 3 |

### Source catalog

| Key | Type | Description | Used in Stage |
|-----|------|-------------|---------------|
| `source_catalog_path` | str | Path to the input source catalog FITS file. Must contain at minimum `ID`, `RA`, `DEC` columns. For Stage 8, also requires `z_grism` and `name_line_exp`. | 5, 6, 7, 8 |

### Data directories and file patterns

| Key | Type | Default | Description | Used in Stage |
|-----|------|---------|-------------|---------------|
| `data_dir` | str | `"."` | Root directory; sub-directories default to paths under this | all |
| `grism_data_dir` | str | `<data_dir>/grism_raw` | Directory containing LW grism `rate.fits` files (Stage-1 output) | 2a |
| `grism_filename_fmt` | str | `"jw*nrc[ab]long_rate.fits"` | Glob pattern to select grism rate files within `grism_data_dir` | 2a |
| `direct_image_dir` | str | `<data_dir>/direct_imaging` | Directory containing NIRCam SW direct-imaging `rate.fits` files | 3, 4 |
| `direct_image_filename_fmt` | str | `"jw*_nrc[a-b][1-4]_rate.fits"` | Glob pattern to select SW rate files within `direct_image_dir` | 3, 4 |
| `calibrated_dir` | str | `<data_dir>/grism_cal` | Output directory for `lv1.5` files, plots, astrometry, and POM catalogs | 2a–5 |
| `extract_dir` | str | `<data_dir>/extract_2d` | Output directory for extracted 2D/1D spectra and diagnostic PDFs | 6, 7, 8 |

### Calibration support data

| Key | Type | Description | Used in Stage |
|-----|------|-------------|---------------|
| `cali_support_dir` | str | Directory containing F. Sun's calibration files: dispersion solutions, sensitivity curves, median super-sky images, POM spectral coverage maps, and astrometry reference tables | 2b, 2c, 5 |
| `default_POM_trans_dir` | str | Directory containing POM transmission maps (`NIRCAM_LW_POM_Mod[AB]_trans.fits`) | 5 |
| `astrometry_cal_table` | str | Path to a pre-computed astrometry correction table. If the file exists it is loaded directly and Stage-4 DAOFIND is skipped; otherwise the table is generated and saved here | 4 |
| `astrometry_ref_table` | str | Path to an external astrometric reference catalog (e.g. Gaia-calibrated HST source list) used for cross-matching in Stage 4 | 4 |

### Direct-image mosaic (for 1D diagnostic plots)

| Key | Type | Description | Used in Stage |
|-----|------|-------------|---------------|
| `image_mosaic_dir` | str | Directory containing large JWST image mosaic FITS files. Used to cut out direct-image thumbnails for the 1D diagnostic PDF. Leave empty to skip mosaic thumbnails. | 7 |
| `image_mosaic_filename_fmt` | str | `str % (field, band)` format string resolving to a mosaic filename inside `image_mosaic_dir`. First placeholder is `image_mosaic_field`, second is the filter name in lower-case. | 7 |
| `image_mosaic_rgb_bands` | list | Three filter names `[blue, green, red]` used to compose the RGB thumbnail in the 1D diagnostic PDF. Example: `[F115W, F200W, F444W]`. | 7 |
| `image_mosaic_field` | str | Field identifier inserted into `image_mosaic_filename_fmt` (e.g. `goods-n`, `goods-s`). | 7 |

### Observation parameters

| Key | Type | Default | Description | Used in Stage |
|-----|------|---------|-------------|---------------|
| `grism_filter` | str | `"F444W"` | NIRCam blocking filter for the grism observations. Determines the wavelength range and which calibration files are loaded. Supported values: `F356W`, `F444W`. | 2a–8 |
| `n_procs` | int | `4` | Number of parallel worker processes for pool-based stages. Set to the number of available CPU cores for best throughput. | 2a–8 |

### Extraction parameters

| Key | Type | Default | Description | Used in Stage |
|-----|------|---------|-------------|---------------|
| `aperture_pix` | float | `15.0` | Cross-dispersion aperture half-width in pixels for 2D spectral extraction. The total aperture is `2 × aperture_pix` pixels wide. Larger values capture more flux from extended sources but increase noise and contamination. | 6 |
| `sigma_hot` | float | `20.0` | Sigma threshold for hot-pixel detection. A pixel is flagged if its value deviates from its neighbours by more than `sigma_hot` times the local noise. Higher values are more conservative (fewer pixels flagged). | 2d |
| `bunit_spec2d` | str | `"DN/s"` | Brightness unit for the co-added 2D spectra written to `extract_dir`. Options: `"DN/s"` (native detector rate units) or `"mJy"` (flux density after applying the instrument sensitivity curve and pixel scale). | 6 |

### Processing control

| Key | Type | Default | Description | Used in Stage |
|-----|------|---------|-------------|---------------|
| `overwrite` | bool | `false` | If `false`, Stages 2a–2c skip files that have already been processed (lv1.5 or background files exist on disk). Set to `true` to force reprocessing. | 2a, 2b, 2c |
| `overwrite_spec2d` | bool | `false` | If `false`, Stage 6 skips sources for which a co-added 2D spectrum file already exists. Set to `true` to re-extract. | 6 |

### Emission-line cutout co-addition (Stage 8)

| Key | Type | Default | Description | Used in Stage |
|-----|------|---------|-------------|---------------|
| `coadd_method` | str | `"simple"` | Co-addition method for emission-line 2D cutouts. `"simple"`: stack in native detector frame, no re-sampling (fast). `"drizzle"`: re-sample each exposure onto a common equatorial grid (produces north-up cutouts at configurable sub-pixel resolution). | 8 |
| `cutout_size_simple` | int | `51` | Side length (pixels) of the cutout in `simple` mode. | 8 |
| `cutout_size_drizzle` | int | `51` | Side length (pixels) of the output grid in `drizzle` mode. | 8 |
| `finalscale_drizzle` | float | `0.035` | Output pixel scale in arcsec/pixel for drizzle co-addition. Smaller values give finer sampling (at the cost of lower S/N per pixel). | 8 |
| `pixfrac_drizzle` | float | `0.8` | Drop-size fraction (pixfrac) passed to the Drizzle algorithm. Values close to 1 produce smoother images; smaller values reduce correlated noise at the cost of more unsampled pixels. | 8 |
| `psf_oversample` | int | `4` | Super-sampling factor for `stpsf` PSF models used as weights in drizzle co-addition. | 8 |

### Minimal working example

```yaml
pid: 9999
grism_filter: F444W
n_procs: 8

source_catalog_path: /data/my_catalog.fits

data_dir: /data/my_survey
grism_data_dir:   /data/my_survey/grism_F444W
direct_image_dir: /data/my_survey/direct_imaging
calibrated_dir:   /data/my_survey/grism_cal
extract_dir:      /data/my_survey/extract

cali_support_dir:      /path/to/nircam_wfss/data/FSun_cal
default_POM_trans_dir: /path/to/nircam_wfss/data/GRISM_NIRCAM
astrometry_cal_table:  /data/my_survey/grism_cal/astrometry.dat
astrometry_ref_table:  /path/to/nircam_wfss/data/FSun_cal/goods_charge_F160W_daofind.cat

image_mosaic_dir: /data/jades/mosaics
image_mosaic_filename_fmt: "hlsp_jades_jwst_nircam_%s_%s_v5.0_drz.fits"
image_mosaic_rgb_bands: [F115W, F200W, F444W]
image_mosaic_field: goods-n

aperture_pix: 15.0
sigma_hot: 20.0
bunit_spec2d: "DN/s"
overwrite: false
overwrite_spec2d: false

coadd_method: "simple"
cutout_size_simple: 31
```

---

## 6. Package Layout

```
nircam_wfss/
├── pyproject.toml              # package metadata and CLI entry point
├── configs/
│   ├── PID1895_FRESCO_GDN_F444W.yaml   # FRESCO GOODS-N F444W example
│   └── PID3577_CONGRESS_GDN_F356W.yaml # CONGRESS GOODS-N F356W example
├── data/
│   ├── FSun_cal/               # dispersion, sensitivity, super-sky, astrometry files
│   ├── GRISM_NIRCAM/           # POM transmission maps
│   └── jwst_pipeline_config/   # JWST pipeline configuration files
└── src/                        # Python source (installed as the nircam_wfss package)
    ├── run_pipeline.py         # CLI entry point — orchestrates all stages
    ├── config.py               # PipelineConfig dataclass and shared constants
    ├── background.py           # Stages 2b & 2c: super-sky and background subtraction
    ├── preprocessing.py        # Stages 2d & 2e: hot-pixel rejection and continuum subtraction
    ├── imaging.py              # Stage 3 & 4: SW image reduction and astrometry calibration
    ├── dispersion.py           # Grism dispersion and spectral trace polynomials
    ├── pom.py                  # Stage 5: POM vignetting and per-frame source catalogs
    ├── extraction.py           # Stages 6, 7, 8: 2D/1D spectral extraction
    ├── noise.py                # 1/f noise subtraction utilities
    ├── wcs_utils.py            # JWST focal-plane ↔ sky coordinate transforms
    └── plotting.py             # Shared matplotlib utilities
```
