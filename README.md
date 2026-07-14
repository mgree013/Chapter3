# Chapter3

Title: Effects of fish on stream macroinvertebrate and lake zooplankton communities

Status of MS: In Prep

Authors: Matthew D. Green, David B. Herbst, and Kurt E. Anderson


# Structure of Code

1) `code/new_lake_analysis.R` - download and process lake zooplankton and environmental data
2) `code/new_lake_analytics.R` - build lake analysis objects and construct figure components
3) `code/new_lake_models.R` - fit lake models for manuscript tables and summaries
4) `code/new_lake_plotting.R` - save explicit manuscript figure PNGs with stable names
5) `code/Plotting.R` - legacy figure export script for the original manuscript workflow

# Data

1) Zooplankton data (Knapp et al. 2020) used in this study can be accessed from Environmental Data Initiative at https://portal.edirepository.org/nis/mapbrowse?scope=edi&identifier=577&revision=2. 

2) Macroinvertebrate data (Green et al. 2022) can be accessed from Dryad at https://doi.org/10.5061/dryad.2fqz612qw.


# Results

### Current generated lake figures
- `fig2a.png` — Lake diversity by fish presence: Shannon diversity
- `fig2b.png` — Lake diversity by fish presence: beta-diversity (LCBD)
- `fig2c.png` — Lake diversity by fish presence: community size
- `fig2d.png` — Lake diversity by fish presence: community biomass
- `fig3a.png` — Lake diversity along elevation: Shannon diversity
- `fig3b.png` — Lake diversity along elevation: beta-diversity (LCBD)
- `fig3c.png` — Lake diversity along elevation: community size
- `fig3d.png` — Lake diversity along elevation: community biomass
- `fig2_panel.png` — combined lake figure panel for Fig 2
- `fig3_panel.png` — combined lake figure panel for Fig 3
- `supp1.png` — supplemental fish presence / elevation figure

### Notes
- The zooplankton density calculation is implemented in `code/new_lake_analysis.R` as `Counts * sample_vol / (zoo_tow_depth * zoo_tow_number * 33.02)`.
- Figure outputs are saved to the project root when `code/new_lake_plotting.R` is run.

# Supplemental

### Fig S1: Fish presence as a function of elevation (m) for lakes and streams.
`code/new_lake_plotting.R` saves `supp1.png` for the current lake supplemental figure.

