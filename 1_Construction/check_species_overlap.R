#check_species_overlap.R
setwd("C:/Users/wande/Documents/EEB498_R/island_ldg")

#dependencies
library(tidyverse)
packageVersion("tidyverse") # ‘2.0.0’
library(stringr)
packageVersion("stringr") # ‘1.5.1’
library(RANN)
packageVersion("RANN") # ‘2.6.2’
library(sf)
packageVersion("sf") # ‘1.0.17’
library(geodist)
packageVersion("geodist") # ‘0.1.0’
library(taxize)
packageVersion("taxize") # ‘0.9.100.1’
library(purrr)
packageVersion("purrr") # ‘1.0.2’

# load in mainland-island pairs and trait data
antdat <- readRDS("data/GABII_islands_efn_ant_trait_data.RDS")
island_mainland_pairs <- readRDS("data/island_mainland_pairs_with_antisland_names.RDS")
antdat.ml <- readRDS("data/GABI_mainlands_efn_trait_data.RDS")

# load in GIFT species checklists and geo information
#Link in the extended geological data to determine oceanic/non-oceanic islands

load("data/GIFT_LatGrad_extended.RData")
species.dat <- checklists
geo.dat <- readRDS("data/geo.dat.extended_foranalyses.rds")

#bring in species dat

species.geo <- species.dat[[1]] %>% as.data.frame() %>%
  select(c('entity_ID', 'entity_class')) %>%
  mutate_at(c('entity_class'), as.character) %>%
  mutate_at(c('entity_ID'), as.integer) 

geo <- geo.dat %>% 
  rename(GDP_satellite = "mean_GDP_satellite", CHELSA_annual_Prec = "mean_wc2.0_bio_30s_12", CHELSA_annual_mean_Temp = "mean_wc2.0_bio_30s_01", 
         Popdensity = "mean_gpw-v4_popdensity_UN_2015", Human_Footprint = "mean_HFP2009_unproj", Human_Influence = "mean_hii_v2geo", 
         elev_range = "range", age_Ma = "age_MA") %>%
  select(c('entity_ID', 'geo_entity', 'area','longitude','latitude', 'biome','dist', 'CHELSA_annual_mean_Temp', 'CHELSA_annual_Prec',
           'GDP_satellite', 'Popdensity', 'Human_Footprint', 'Human_Influence', 'elev_range', 'geology', 'age_Ma')) %>%
  mutate_at(c('area','longitude','latitude', 'dist','CHELSA_annual_mean_Temp', 'CHELSA_annual_Prec',
              'GDP_satellite', 'Popdensity','Human_Footprint', 'Human_Influence', 'elev_range', 'age_Ma'), as.numeric) %>%
  mutate_at(c('geo_entity','biome','geology'), as.character) %>%
  mutate_at(c('entity_ID'), as.integer) %>%
  left_join(species.geo, by = "entity_ID") %>%
  distinct(entity_ID, .keep_all = TRUE) %>%
  mutate(entity_class = case_when(entity_class == "Island" ~ "Island",
                                  entity_class == "Island Part" ~ "Island",          
                                  entity_class == "Island Group" ~ "Island",
                                  entity_class == "Mainland" ~ "Mainland",
                                  entity_class == "Island/Mainland" ~ "undetermined")) %>%
  mutate(geology = case_when(geology == "atoll" ~ "dev", 
                             geology == "floor" ~ "dev", 
                             geology == "floor/volcanic" ~ "dev",
                             geology == "volcanic" ~ "dev",
                             geology == "shelf" ~ "nondev",
                             geology == "fragment" ~ "nondev",
                             geology == "atoll/shelf/volcanic" ~ "nondev",
                             geology == "fragment/shelf/volcanic" ~ "nondev")) %>% 
  mutate(entity_class2 = case_when(geology == "dev" ~ "Oceanic",                      
                                   geology == "nondev" ~ "Non-oceanic",
                                   entity_class =="Mainland" ~ "Mainland"))

### bind trait data to plant species and get list by geo
#Plant Trait Data: All Regions
#Binding EFN presence/absence trait data to GIFT database
#bring in EFN species data
#to construct data/efn_resolved_full.csv, run 1_Construction/EFNs_taxize.Rmd

efndat <- read.csv("data/efn_resolved_full.csv", header = TRUE) %>% 
  rename(species = matched_name) %>% 
  mutate("species" = ifelse(is.na(species), Phy, species)) %>% 
  select("EFN", "species") %>% 
  distinct(species, .keep_all = TRUE)

#with species checklist data, keep native, assign EFN using species, fill 0 if no info is available
#join with geo data and keep unique sp by loc (entity_ID)

occurrences <- readRDS("data/species_export_GBIF_filtered.RDS")

occurrences <- occurrences %>% 
  select(c("species", "occ_count"))

species <- species.dat[[2]] %>% as.data.frame() %>% 
  mutate(species = paste(genus, species_epithet, sep = " ")) %>% 
  select(c("entity_ID", "family", "native", "naturalized", "species", "name_ID", "genus")) %>%
  mutate_at(c('family', 'species', 'name_ID', 'genus'), as.character) %>%
  mutate_at(c('native', 'naturalized'), as.numeric) %>%
  mutate_at(c('entity_ID'), as.integer) %>%
  select(-c("native", "naturalized")) %>%
  left_join(efndat, by = c("species")) %>% 
  mutate(EFN = ifelse(is.na(EFN), 0, EFN)) %>% 
  left_join(occurrences, by = c("species")) %>% 
  filter(occ_count >= 100) %>% 
  drop_na(EFN) %>%
  left_join(geo, by = "entity_ID") %>%                                          
  group_by(entity_ID) %>% 
  distinct(species, .keep_all = TRUE) %>%             
  ungroup() 

# Plant trait data: require intersection between mainland and island species to count island species presence
islandIDs <- unique(island_mainland_pairs$entity_ID.i)

# Step 1: for each row of island_mainland_pairs, extract matching species
sp_common_all.plant <- island_mainland_pairs %>%
  rowwise() %>%
  mutate(
    shared_species = list(
      {
        island_sp <- species %>%
          filter(entity_ID == entity_ID.i) %>%
          pull(species)
        
        mainland_sp <- species %>%
          filter(entity_ID == entity_ID.ml) %>%
          pull(species)
        
        intersect(island_sp, mainland_sp)
      }
    )
  ) %>%
  unnest(shared_species) %>%
  rename(species = shared_species) %>%
  mutate(
    entity_ID = entity_ID.i,
    sp_on_ml = "Y"
  ) %>%
  select(pairID, entity_ID, species, sp_on_ml)

sp.overlap.stats <- sp_common_all.plant %>% 
  group_by(pairID) %>% 
  summarize(total = n())

summary(sp.overlap.stats$total)

species.overlap <- species %>% 
  left_join(sp_common_all.plant, by = c("entity_ID", "species")) %>% 
  filter(
    (entity_ID %in% islandIDs & sp_on_ml == "Y") |
    !(entity_ID %in% islandIDs)
  )

sprich.o.efn.pres <- species %>%
  filter(EFN == 1) %>%
  group_by(entity_ID) %>%                                                                                           
  summarise(efn.pres = n()) 

sprich.o.efn.abs <- species %>%
  filter(EFN == 0) %>%
  group_by(entity_ID) %>%                                                                                           
  summarise(efn.abs = n()) 

sprich.efn.pres <- species %>%
  filter(EFN == 1) %>%
  group_by(entity_ID) %>%                                                                                           
  summarise(efn.pres = n()) 

sprich.efn.abs <- species %>%
  filter(EFN == 0) %>%
  group_by(entity_ID) %>%                                                                                           
  summarise(efn.abs = n()) 