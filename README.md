mammal-abundance
====================

[![](https://zenodo.org/badge/DOI/10.5281/zenodo.154200.svg)](http://doi.org/10.5281/zenodo.154200)

An N-mixture model approach for estimating abundance from encounter histories of untagged animals.

Specifically, this analysis examined the short-term effects of silvicultural treatments on abundance of common small mammal species at the Hardwood Ecosystem Experiment (HEE) sites. Trapping data were collected yearly from 32 permanent grids established in 2007, yielding 2 pre- and 3 post-harvest years.

There are several sub-analyses included in this repository:

1. Small mammal abundance - abundance estimates from an N-mixture model for four small mammal species (white-footed mouse, eastern chipmunk, short-tailed shrew and pine vole) based on repeated count data.

2. A validation of the N-mixture approach, by comparing abundance estimates with mark-release-recapture (MRR) estimates at the same sites, and via simulation.

3. A spatially explicit capture-recapture (SCR) analysis conducted at a subset of trapping grids were animals were marked, with the goal of estimating the effective sample area of the grids.

4. Brief analysis of microhabitat data collected at each trap location.

The results (pre vs. post-harvest small mammal abundance) are published in the following article:

[*Kellner, Kenneth F.; Urban, Natasha A.; Swihart, Robert K. 2013. Short‐term responses of small mammals to timber harvest in the United States Central Hardwood Forest Region. Journal of Wildlife Management 77:1650-1663*](http://onlinelibrary.wiley.com/doi/10.1002/jwmg.613/full)

Metadata
====================

The *analysis_[ ].R* and *sim_abundance.R* files contain the framework of the Bayesian analysis in JAGS for each of the project components listed above. Additional information is provided in the comments for each file.

The *models* subfolder contains the model description (in the BUGS language) for each JAGS analysis contains the corresponding JAGS analysis outputs in Rdata objects. The *figures* folder contains code generating the figures in the final manuscript.

The *data* subfolder contains example CSV files representative of the actual data used in this study. I do not currently have permission to publish the complete raw data. Columns in each CSV file are described below.

1.  *example_captures.csv* contains the raw small mammal capture records from a single year and have been processed to fix errors. Each row represents a single capture/trap record and the following columns are available:  
    -  Grid - the trapping grid ID number
    -  Check - the check day/time identification (5 days with 2 checks/day; 1P = day 1, PM; 3A = day 3, AM, etc.)
    -  Trap - the trap ID; each grid had 3 rows designated by letters with 9 traps in each row designated by numbers.
    -  Species - the small mammal species captured, if applicable. Represented by a code (WFMO = white-footed mouse, EACH = eastern chipmunk, SHSH = short-tailed shrew, PIVO = pine vole).
    -  id - the tag ID for the animal, if it was eartagged.
    -  Fate - The fate of the animal and/or trap, most impmortantly '1' for new capture, '2' for re-capture, and '7' or '8' for missing/disturbed trap, etc.

2.  *example_grid_covs.csv* contains information about each trapping grid location:
    -  temp1-5 - Mean daily temperature for each of the 5 trapping days (measured on-site).
    -  jd1-5 - Julian day for each of the trapping days.
    -  aspect - The site aspect; 1 for northeast and 0 for southwest.
    -  mast1-5 - Standardized index of mast availability at the site based on nearby acorn collection data.
    -  precl-post2 - Treatment identification indices.

3.  *example_microsite_means.csv* contains mean values for microhabitat variables at each trapping grid location (summarized over all trap microsite locations):
    -  grid - the trapping grid ID number.
    -  herb - percent herbaceous plant cover.
    -  wood - percent woody plant cover.
    -  cwd - length of coarse woody debris in a 1m radius around trap location (m).
    -  litter - litter depth (cm).
    -  year - study year (1-5).
    -  treat - silvicultural treatment code.


