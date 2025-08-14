## Overview

This package contains the R files and data necessary to  create the species-level and island-level analyses.

### References
Delavaux, C. S., T. W. Crowther, J. D. Bever, P. Weigelt, and E. M. Gora. 2024. Mutualisms weaken the latitudinal diversity gradient among oceanic islands. Nature (London) 627:335-339. https://doi.org/10.1038/s41586-024-07110-y
https://github.com/c383d893/IslandLDG

Pichler, M., & Hartig, F. (2024). Is there a robust effect of mainland mutualism rates on species richness of oceanic islands? arXiv preprint arXiv:2411.15105. 
https://github.com/MaximilianPi/IslandLDG

König, C., Weigelt, P., Taylor, A., Stein, A., Dawson, W., Essl, F., Pergl, J., Pyšek, P., Van Kleunen, M., & Winter, M. (2021). Source pools and disharmony of the world's island floras. Ecography, 44(1), 44-55. 
[https://doi.org/10.1111/ecog.05174](https://doi.org/10.1111/ecog.05174)
https://doi.org/10.5061/dryad.6wwpzgmwm

## Data Sources
#### GIFT data:  GIFT R package
Retrieved and processed via code provided in 
Delavaux, C. S., T. W. Crowther, J. D. Bever, P. Weigelt, and E. M. Gora. 2024. Mutualisms weaken the latitudinal diversity gradient among oceanic islands. Nature (London) 627:335-339. https://doi.org/10.1038/s41586-024-07110-y
https://github.com/c383d893/IslandLDG  

files:  

- data/GIFT_LatGrad_extended.RData (build with 1_Construction/GIFT_Extraction.R) - Geographic regions including islands and mainlands, floral checklists
- data/GIFT_LatGrad_extended_foranalyses.RData (sourced directly from https://github.com/c383d893/IslandLDG) - geology information for determining oceanic vs non-oceanic islands`
- data/GIFT_extra_traits_extended.Rds - woodiness, lifecycle, growth form

#### GABI
Guenard, B., M. Weiser, K. Gomez, N. Narula, and E. P. Economo. 2017. The Global Ant Biodiversity Informatics (GABI) database: synthesizing data on the geographic distribution of ant species (Hymenoptera: Formicidae). Myrmecological News 24:83-89.

Referenced in:
Nathan, P., E. P. Economo, B. Guénard, A. K. Simonsen, and M. E. Frederickson. 2023. Generalized mutualisms promote range expansion in both plant and ant partners. Proceedings of the Royal Society B: Biological Sciences 290:20231083.
Available at: https://doi.org/10.5061/dryad.jm63xsjh6

Contains: Release package of the GABI database plus shapefiles corresponding to the geographic regions referenced in GABI records

#### GABI-I
Liu, C., E. P. Economo, and B. Guénard. 2023. GABI‐I: The global ant biodiversity informatics‐island database. Wiley Online Library.
Download link directly embedded in paper.  

#### EFN Plant Trait Data
Weber, M., L. Porturas, and K. Keeler. 2015. World list of plants with extrafloral nectaries.
http://www.extrafloralnectaries.org/  
- data/efn_full_database.csv

Processed according to / sourced from:
Nathan, P., E. P. Economo, B. Guénard, A. K. Simonsen, and M. E. Frederickson. 2023. Generalized mutualisms promote range expansion in both plant and ant partners. Proceedings of the Royal Society B: Biological Sciences 290:20231083.
https://github.com/The-Frederickson-Lab/ant-legume-range  

files: 
- data/efn_resolved.csv

#### Ant EFN visiting Data
Kaur, K. M., P.-J. G. Malé, E. Spence, C. Gomez, and M. E. Frederickson. 2019. Using text-mined trait data to test for cooperate-and-radiate co-evolution between ants and plants. PLoS Computational Biology 15:e1007323.
https://doi.org/10.5061/dryad.3630t8j  

files: 
- data/Species_EFN_Data.csv

#### Island-Mainland Pairings
König, C., Weigelt, P., Taylor, A., Stein, A., Dawson, W., Essl, F., Pergl, J., Pyšek, P., Van Kleunen, M., & Winter, M. (2021). Source pools and disharmony of the world's island floras. Ecography, 44(1), 44-55. 
[https://doi.org/10.1111/ecog.05174](https://doi.org/10.1111/ecog.05174)
files: 
- data/konig - populate from https://doi.org/10.5061/dryad.6wwpzgmwm
- 1_Construction/source_pool_estimation.R
- 1_Construction/source_pool_estimation_ants.R

#### Ant Phylogeny
Nelsen, M. P., Ree, R. H., & Moreau, C. S. (2018). Ant–plant interactions evolved through increasing interdependence. Proceedings of the National Academy of Sciences, 115(48), 12253-12258. 
https://datadryad.org/dataset/doi:10.5061/dryad.ft4sn88  

files: 
- data/ant_tree.tre


## R Files

### Constructing the main analysis (Island-Mainland pairs)

#### Data preparation

To construct analysis files, proceed in this order:
1_Construction/GIFT_Extraction.R - initial API request for Tracheophytes from GIFT - note this file takes about an hour to run, produces a large RData file data/GIFT_LatGrad_extended.RData

1_Construction/GIFT_Restricted_For_Extended.R - this interfaces with the restricted API endpoint of GIFT, to request geological information.  The output of this file has been directly provided in the code repo as data/geo.dat.extended_foranalyses.rds

1_Construction/EFNs_taxize.Rmd - run to produce data/efn_resolved.csv
input: data/efn_full_database.csv (provided)

1_Construction/GBIF_occurrence_data.Rmd - run to produce data/species_export_GBIF_filtered.RDS, data/ant_i_species_export_GBIF_filtered.RDS, data/ant_species_export_GBIF_filtered.RDS

1_Construction/source_pool_estimation.R - run to produce data/konig/final_predicted_mainland_island_pairs.Rds

1_Construction/source_pool_estimation_ants.R - run to produce data/konig/pairs_pred_ants.RData

1_Construction/Island_Mainland_and_Trait_Construction.Rmd
This file has several purposes
- to geographically cross-correlate islands between GABI-I and GIFT datasets, where islands have slightly different reported latitude/longitude coordinates
- to match each unique island in GIFT with its most likely source mainland following Konig et al (2022)
- to add EFN presence/absence data to GIFT floral checklists
- to add EFN visitor/nonvisitor data to GABI and GABI-I ant species occurrence records
- to filter species checklists by GBIF occurrence counts
- to create a many-to-many index that binds island name from GABI with entity_ID of an island in GIFT and entity_id of a nearest mainland in GIFT

inputs:   
- data/GABI-I.csv as downloaded above
- data/GIFT_LatGrad_extended.RData as constructed above
- data/geo.dat.extended_foranalyses.rds as provided in repo
- data/Species_EFN_Data.csv as downloaded above
- data/efn_resolved_full.csv as constructed above
- GABI_shapefiles/GABI_Data_Release1.0_18012020.csv as downloaded above
- GABI_shapefiles\\Bentity2_shapefile_fullres.shp as downloaded above
- data/species_export_GBIF_filtered.RDS as constructed above
- data/ant_i_species_export_GBIF_filtered.RDS as constructed above
- data/ant_species_export_GBIF_filtered.RDS as constructed above
- data/konig/final_predicted_mainland_island_pairs.Rds as constructed above
- data/konig/pairs_pred_ants.RData as constructed above

outputs:   
- data/island_mainland_pairs_with_antisland_names.RDS - the many to many index
- data/island_mainland_pairs_ants.RDS 
- data/GABIGABI_efn_trait_data.RDS
- data/GIFT_EFN_trait_data.RDS

#### Analysis
run file 2_Analyses\Island_Mainland_Pairs_Analysis.Rmd after running through construction steps above

### Requiring Species Overlap
To recreate the analysis which requires species overlap between islands and source mainlands, use 
1_Construction\Island_Mainland_and_Trait_Construction-RequireSpeciesMatch.Rmd

Then, run 2_Analysis\Island_Mainland_Pairs_Analysis.Rmd with the boolean at the top set to TRUE.

### PGLS Model - Plants
Data source:  data/efn_resolved_full.csv, data/GIFT_LatGrad_extended.RData, data/geo.dat.extended_foranalyses.rds, data/GIFT_extra_traits_extended.Rds as constructed above

Construction/analysis: 1_Construction/Construction_PGLS.Rmd
Note: this is a very large phylogeny and takes a few hours to run

### PGLS Model - Ants
Data source: data/GABI-I.csv, data/Species_EFN_Data.csv, GABI_shapefiles/GABI_Data_Release1.0_18012020.csv, data/ant_tree.tre

Construction/analysis: 1_Construction/Construction_PGLS_ants.Rmd