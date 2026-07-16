# Chapter3

Title: Trout have weak effects on zooplankton diversity but strong effects on community biomass

Status of MS: Limnology and Oceanography 

Authors: Matthew D. Green and Kurt E. Anderson

Link: https://aslopubs.onlinelibrary.wiley.com/doi/10.1002/lno.70041

# How to Reproduce

Run the scripts **in order** from the project root:

```r
setwd("/Users/matthewgreen/Library/CloudStorage/Dropbox/Manuscipts/Chapter 3/Chapter3")
source("code/new_lake_analysis.R")   # ~5–10 min (downloads EDI data)
source("code/new_lake_analytics.R")  # builds all figures + models + saves PNGs
```

`new_lake_plotting.R` is now a thin wrapper that can optionally re-save figures from a loaded RData, but all figures are already saved by `new_lake_analytics.R`.


# Structure of Code

| Script | Purpose |
|--------|---------|
| `code/Download_Slip_Data.R` | Downloads all 13 EDI data tables (edi.577.2) |
| `code/new_lake_analysis.R` | Processes raw data → `data/new_lake_processed.RData` |
| `code/new_lake_analytics.R` | Builds all manuscript figures + models → `data/new_lake_analytics.RData` + saves PNGs to `Newfigs/` |
| `code/new_lake_plotting.R` | Optional: re-save figures from `new_lake_analytics.RData` |
| `code/Lake.analysis.R` | Full legacy + current analysis (used interactively) |
| `code/Plotting.R` | Legacy figure export script |


# Density Calculation

Zooplankton density (individuals/L):

```
density = (Counts / Max.Subsample) × sample_vol
          ─────────────────────────────────────────
          zoo_tow_depth × zoo_tow_number × 33.02
```

- `sample_vol` comes from `dt13` (lab processing table) joined via **`collect_date`** (not `survey_date`, which is NA for most dt13 rows).
- Diversity metrics (N0, Shannon H, N1, LCBD) are computed from **subsample-corrected counts** so that the `sample_vol` join never causes NA diversity even for the ~128 sites where dt13 has no matching record.
- `Com.Size` (total density) uses the density-based matrix and will be 0 for those 128 sites.


# Data

Zooplankton data (Knapp et al. 2020): <https://portal.edirepository.org/nis/mapbrowse?scope=edi&identifier=577&revision=2>

# Output Figures (saved to `Newfigs/`)

| File | Manuscript Figure | Description |
|------|-------------------|-------------|
| `Figure1.png` | Fig 1 (a–d) | Shannon diversity + LCBD × fish presence (boxplots) and elevation (scatter) |
| `Figure2.png` | Fig 2 | Beta-diversity components: β_bal (turnover), β_gra (nestedness), β_bray (total) |
| `Figure3.png` | Fig 3 | Relative change in species density vs log body mass |
| `Figure4.png` | Fig 4 | Community weighted mean (CWM) body mass by fish presence |
| `FigureS1.png` | Fig S1 | Fish presence vs elevation + lake site counts |
| `FigureS2.png` | Fig S2 | Species log density (a) + proportion of lakes occupied (b), ordered by body size |
| `FigureS3.png` | Fig S3 | Body mass of species absent from fish or fishless sites |
| `FigureS4.png` | Fig S4 | CWM vs elevation, faceted by fish presence |


# Models (saved in `data/new_lake_analytics.RData`)

| Object | Table | Response | Predictors |
|--------|-------|----------|------------|
| `m_N1_*` + `table1_N1` | Table 1 | Shannon diversity (N1) | elevation, fish, interaction, null |
| `m_LCBD_*` + `table1_LCBD` | Table 1 | LCBD | elevation, fish, interaction, null |
| `m_beta_*` + `table2_beta` | Table 2 | β diversity component | bal vs gra |
| `m_density_*` + `table3_density` | Table 3 | Relative density change | log body mass, null |
| `m_CWM_*` + `table4_CWM` | Table 4 | CWM | fish × elevation, elevation, fish, null |
| `m_absent_*` + `tableS3_absent` | Table S3 | Log body mass | absent from fish vs fishless |


