#Based on Konig et al. (2020) in Ecography
#See https://datadryad.org/dataset/doi:10.5061/dryad.6wwpzgmwm

library(tidyverse)
library(vegan)
library(gdm)
library(geosphere)
library(geodist)
library(raster)
library(cowplot)
library("rnaturalearth")
library("rnaturalearthdata")
library(ggspatial)
library(segmented)
library(sf)

#Load meta-data for 178 islands and 735 mainlands
load("data/konig/geoentities_meta.RData")
#Subset metadata to islands only
meta_islands <- geoentities_meta[geoentities_meta$entity_class == "Island", ]

#Load GABI-I data
gabi_i <- read.csv("data/GABI-I.csv", header = TRUE)
gabi_i$species <- paste0(gabi_i$Genus.Name.Pub, " ", gabi_i$Species.Name.Pub) #Species name

# Function to calculate UTM zone based on longitude
get_utm_zone <- function(longitude) {
  zone <- floor((longitude + 180) / 6) + 1
  return(zone)
}

#Make a list of ant islands
antislands <- gabi_i %>% 
  mutate(latlong = paste(Lat, Long, sep = ",")) %>% 
  distinct(latlong, .keep_all = TRUE) %>%
  dplyr::select(c("Island.Name", "Bentity2.Name", "latlong")) %>% 
  separate(col = latlong, into = c("lat", "long"), sep=",")

antislands$lat <- as.numeric(antislands$lat)
antislands$long <- as.numeric(antislands$long)

# set UTM zone for each lat/long combination
antislands$utm_zone <- sapply(antislands$long, get_utm_zone)

# handle northern vs. southern hemisphere
sf_points <- st_as_sf(antislands, coords = c("long", "lat"), crs = 4326)

# Check for NA or invalid values in the coordinate columns
sum(is.na(meta_islands$latitude) | is.na(meta_islands$longitude))  # Count NA values in islands
sum(is.na(antislands$lat) | is.na(antislands$long))  # Count NA values in antislands

# Remove rows with NA or invalid values in both data frames
islands_clean <- meta_islands[!is.na(meta_islands$latitude) & !is.na(meta_islands$longitude), ]
antislands_clean <- antislands[!is.na(antislands$lat) & !is.na(antislands$long), ]

# Calculate the distances between all pairs of points
dist_matrix <- geodist(
  x = antislands_clean[, c("long", "lat")],
  y = islands_clean[, c("longitude", "latitude")],
  measure = "haversine" # Great-circle distance
)

# Find the nearest neighbors (minimum distance for each query point)
nn_indices <- apply(dist_matrix, 1, which.min)
nn_distances <- apply(dist_matrix, 1, min)

#Save nearest neighbors and distances
antislands_clean$nearest_neighbour_entity_ID <- islands_clean$entity_ID[nn_indices]
antislands_clean$distance <- nn_distances

#Join with island data
antislands_clean <- antislands_clean %>% 
  left_join(islands_clean, by=c("nearest_neighbour_entity_ID" = "entity_ID"))

#Filter island pairs for distances under 1 km (1000m)
island_pairs <- antislands_clean %>% 
  filter(distance < 1000) %>% 
  mutate(nn_final = nearest_neighbour_entity_ID) # only 142 islands paired based on distance

#Check for duplicates
island_pairs_duplicates <- island_pairs %>% 
  group_by(nearest_neighbour_entity_ID) %>% 
  filter(n() > 1) %>% 
  ungroup()

length(unique(island_pairs_duplicates$geo_entity))

# decide duplicates by shortest distance
island_pairs <- island_pairs %>% 
  group_by(nearest_neighbour_entity_ID) %>% 
  mutate(
    nn_final = case_when(
      n() > 1 & distance == min(distance) ~ nn_final,
      n() > 1 & distance != min(distance) ~ NA,
      n() <= 1 ~ nn_final
    )
  ) %>% 
  ungroup() 
island_pairs <- subset(island_pairs, !is.na(nn_final))

# drop columns for lat and long from GABI-I, use those from GIFT
island_pairs <- island_pairs %>% 
  dplyr::select(-c(lat, long)) %>% 
  rename(entity_ID.i = nearest_neighbour_entity_ID) %>% 
  rename(lat = latitude) %>% 
  rename(long = longitude)

antisland_entity_ids <- unique(island_pairs$entity_ID.i)

#Extract species lists for each island that has been matched between GABI and GIFT
gabi_i_matched <- merge(gabi_i, island_pairs, by="Island.Name")

#Load GABI mainlands data
gabi_mainland <-  read_csv("GABI_shapefiles/GABI_Data_Release1.0_18012020.csv")
gabi_mainland$species <- paste0(gabi_mainland$genus_name_pub, " ", gabi_mainland$species_name_pub) #Species name

#Merge all of GABI
gabi_i_matched$bentity2_name <- gabi_i_matched$Island.Name
gabi_mainland$entity_class <- "Mainland"
gabi_all <- rbind(gabi_i_matched[, c("bentity2_name", "entity_class", "species")], gabi_mainland[, c("bentity2_name", "entity_class", "species")])

#Remove duplicates
gabi_all <- gabi_all[!duplicated(gabi_all), ]

############################################
### Compare GABI polygons and grid cells ###
############################################

#6495 grid cells from Konig et al. 2021
# load global equal area grid and prepared environmental variables for each grid cell
load(file = "data/konig/grid.RData")
load(file = "data/konig/env_grid.RData")

#Load GABI polygons (shapefiles)
Bentity2 <- st_read("GABI_shapefiles/Bentity2_shapefile_fullres.shp")
# Ensure the CRS is WGS84 (latitude/longitude)
Bentity2 <- st_transform(Bentity2, crs = 4326)
# Check geometry validity
validity <- st_is_valid(Bentity2)
# Inspect invalid geometries
invalid_geometries <- Bentity2[!validity, ]
print(invalid_geometries)
# Fix invalid geometries
Bentity2 <- st_make_valid(Bentity2)
# Check validity again
validity <- st_is_valid(Bentity2)
if (all(validity)) {
  print("All geometries are now valid.")
} else {
  print("Some geometries are still invalid.")
}

#World map base layer
world <- ne_countries(scale = "medium", returnclass = "sf")

#Give everything the same coordinate system
polygons_sf <- st_set_crs(st_as_sf(st_make_valid(Bentity2)), "+proj=longlat +datum=WGS84")
polygons_sf <- st_make_valid(polygons_sf)
grid_cells <- st_set_crs(st_as_sf(grid), "+proj=longlat +datum=WGS84")
grid_cells <- st_make_valid(grid_cells)

#Visualize overlap between 6495 grid cells in Konig and GABI polygons
p1 <- ggplot()+
  geom_sf(data=world) +
  xlab("Longitude") + ylab("Latitude") +
  geom_sf(data=polygons_sf, color="black")+
  geom_sf(data=grid_cells, fill="pink", alpha=0.1)+
  theme_cowplot()+
  theme(legend.position="none")
p1 

save_plot("grid_cells_overlayed_on_GABI.pdf", p1)

#Map grid cells to GABI polygons for the purposes of adding environmental data from Konig et al. to each GABI polygon
sf_use_s2(FALSE)
sf_cent <- st_centroid(polygons_sf)
grid_gabi <- st_join(sf_cent, grid_cells, left=TRUE)
env_grid$ID <- gsub("cell_", "", env_grid$entity_ID)
grid_gabi_env <- merge(grid_gabi, env_grid, by.x = "ID", by.y = "ID")
grid_gabi_env$entity_class <- "Mainland"
grid_gabi_env <- grid_gabi_env[c("ID", "BENTITY2_N" , "entity_class", "latitude", "longitude", "area", "T_mean", "T_var", "P_mean", "P_var")]
meta_islands <- meta_islands[c("entity_ID", "geo_entity", "entity_class", "latitude", "longitude", "area", "T_mean", "T_var", "P_mean", "P_var")]
colnames(grid_gabi_env) <- c("entity_ID", "geo_entity", "entity_class", "latitude", "longitude", "area", "T_mean", "T_var", "P_mean", "P_var", "geometry")
meta_islands <- as.data.frame(meta_islands)
meta_mainlands <- as.data.frame(grid_gabi_env) 
meta_mainlands <- meta_mainlands[ , -11]
meta_islands[which(meta_islands$entity_ID %in% meta_mainlands$entity_ID),]
meta_mainlands$entity_ID <- meta_mainlands$entity_ID+20000
geoentities_meta_ants <- rbind(meta_mainlands, meta_islands)

###########################
### Calculate distances ###
###########################
# 1. Compositional dissimilarity (Species level)
gabi_all <- gabi_all[gabi_all$bentity2_name %in% geoentities_meta_ants$geo_entity, ]
gabi_all <- merge(gabi_all, geoentities_meta_ants[, c("entity_ID", "geo_entity")], by.x="bentity2_name", by.y="geo_entity")
spec_tab_ants = unclass(table(gabi_all[,c("entity_ID", "species")])) # look at species composition
#spec_tab is a matrix of sites (422 entities in the GABI-I dataset) by species (39039 species)
spec_dist_ants = betadiver(spec_tab_ants, "sim") #   Turnover / Simpsons Beta diversity
spec_dist_ants = as.matrix(spec_dist_ants) # Convert to matrix

# 2. Geographical distance
geo_dist_ants = round(distm(as.matrix(geoentities_meta_ants[, c("longitude", "latitude")]), fun = distHaversine)/1000,1)+0.1
colnames(geo_dist_ants) = geoentities_meta_ants$entity_ID
rownames(geo_dist_ants) = geoentities_meta_ants$entity_ID

save(list = c("geo_dist_ants", "spec_dist_ants"), file = "data/konig/distances_ants.RData")
gc() #Garbage collection

############################################
###  Generalized dissimilarity modeling  ###
############################################
# We fit a model to reflect the default behaviour of floristic similarity on the mainland 
# using distance + climatic variables. Given new data, we can predict the expected number of shared species
# relative to a global equal area grid in order to localize potential source pools for islands. 
spec_dist_ants = cbind("entity_ID" = as.numeric(rownames(spec_dist_ants)), spec_dist_ants)
geo_dist_ants = cbind("entity_ID" = as.numeric(rownames(geo_dist_ants)), geo_dist_ants)

#This just makes a list of mainland entity IDs and a list of island entity ID
mainlands = unique(gabi_all$entity_ID[gabi_all$entity_class == "Mainland"])
islands = unique(geoentities_meta_ants$entity_ID[geoentities_meta_ants$entity_class == "Island"])

### 1. Model fitting
spec_dist_ml = spec_dist_ants[paste(mainlands), c("entity_ID", paste(mainlands))]
geo_dist_ml = geo_dist_ants[paste(mainlands), c("entity_ID", paste(mainlands))]
env_data_ml = geoentities_meta_ants %>% filter(entity_ID %in% spec_dist_ml[,"entity_ID"]) %>% dplyr::select(entity_ID, longitude, latitude, T_mean, T_var, P_mean, P_var)

#This fits a model using mainland data
gdm_tab_ml = formatsitepair(bioData = spec_dist_ml, bioFormat = 3, siteColumn = "entity_ID", abundance = F, XColumn = "longitude", YColumn = "latitude", predData = env_data_ml, distPreds = list("Geo_Distance" = geo_dist_ml))
gdm_fit_ml = gdm(gdm_tab_ml, geo = F)
summary(gdm_fit_ml)

#Save the model
save(gdm_fit_ml, file = "data/konig/gdm_fit_ml_ants.RData")

### 2. Prepare prediction Data
# Combine environmental data of grid cells and islands 
env_data_is = geoentities_meta_ants %>% filter(entity_ID %in% islands) %>% dplyr::select(entity_ID, longitude, latitude, T_mean, T_var, P_mean, P_var)
#Make entity ID numeric, as it no longer can be a string
env_grid$entity_ID <- gsub("cell_", "", env_grid$entity_ID)

# To distinguish 'real' entities from grid cell entities, add 100000 to grid cell entity IDs
#Any value over 100000 is from the grid
env_grid$entity_ID <- as.numeric(env_grid$entity_ID)+100000
env_grid = rbind(env_data_is, env_grid[, -8])

#Double check that entity ID is numeric
env_grid$entity_ID <- as.numeric(env_grid$entity_ID)

# Prepare geographic distance matrix for gdm function 
geo_dist_pred <- sapply(1:nrow(env_grid), function(i){
  distGeo(as.matrix(env_grid[i ,c("longitude","latitude")]), env_grid[,c("longitude","latitude")])
})/1000
geo_dist_pred = cbind(env_grid$entity_ID, geo_dist_pred)
colnames(geo_dist_pred) = c("entity_ID", env_grid$entity_ID)

# Prepare empty compositional dissimilarity matrix (gdm::formatsitepair requires this)
spec_dist_pred = cbind(env_grid$entity_ID, matrix(0, nrow = nrow(env_grid), ncol = nrow(env_grid)))
colnames(spec_dist_pred) = c("entity_ID", env_grid$entity_ID)

# Prepare gdm table
gdm_tab_pred = formatsitepair(bioData = spec_dist_pred, bioFormat = 3, siteColumn = "entity_ID", abundance = F,
                              XColumn = "longitude", YColumn = "latitude", predData = env_grid, 
                              distPreds = list("Geo_Distance" = geo_dist_pred))
save(gdm_tab_pred, file = "data/konig/gdm_tab_pred_ants.RData")

### 3. Predict
#Changed from predict.gdm to predict, but seems to work fine
gdm_pred = predict(gdm_fit_ml, gdm_tab_pred, geo = F, time = F)
gdm_pred_mat = matrix(NA, nrow = nrow(env_grid), ncol = nrow(env_grid))
gdm_pred_mat[which(lower.tri(gdm_pred_mat))] = gdm_pred
gdm_pred_mat = as.matrix(as.dist(gdm_pred_mat))

entity_names = sort(env_grid$entity_ID)
colnames(gdm_pred_mat) = entity_names
rownames(gdm_pred_mat) = entity_names

save(gdm_pred_mat, file = "data/konig/gdm_pred_mat_ants.RData")

#############################
### Most likely lat/long ###
#############ä###############

#From this prediction matrix, extract the most likely source grid cell for each island

#Make matrix a dataframe
model_pred_df <- as.data.frame(gdm_pred_mat)

#Keep islands only as columns
model_pred_df <- model_pred_df[, paste(islands)]

#Keep mainland grid cells as rows 
model_pred_df <- model_pred_df[which(as.numeric(rownames(model_pred_df)) > 100000), ]

#For each island, find the top source mainlands
pred_latlong <- data.frame(island = NULL, mainland_cell = NULL, turnover = NULL, source_longitude = NULL, source_latitude = NULL, T_mean = NULL, T_var = NULL, P_mean = NULL, P_var = NULL)

for(island in islands){
  tmp <- data.frame(pred = model_pred_df[, paste(island)], grid=rownames(model_pred_df))
  tmp <- merge(tmp, env_grid, by.x="grid", by.y="entity_ID", all.x=T, all.y=F)
  tmp$island <- island
  tmp <- tmp[which(tmp$pred == min(tmp$pred)), ]
  colnames(tmp) <- c("mainland_cell", "turnover", "source_longitude", "source_latitude", "T_mean", "T_var", "P_mean", "P_var", "island")  
  tmp <- tmp[, c(9, 1:8)]
  pred_latlong <- rbind(pred_latlong, tmp)
}

#For any island, plot model predictions
#plot <- ggplot()+
#  geom_sf(data=world) +
#  xlab("Longitude") + ylab("Latitude") +
#  geom_point(data=tmp, aes(y=latitude, x=longitude, color=pred))+
 # theme_cowplot()
#plot

#Add island lat long
island_latlong <- env_grid[which(islands %in% env_grid$entity_ID), c("entity_ID", "longitude", "latitude")]
colnames(island_latlong) <- c("island", "island_longitude", "island_latitude")
pred_latlong <- merge(island_latlong, pred_latlong, by="island")

#A little code to prevent lines from crossing between the South Pacific on the left and Australia on the right
pred_latlong$long_dist <- pred_latlong$source_longitude-pred_latlong$island_longitude
pred_latlong$direction <- ifelse(abs(pred_latlong$long_dist)>180 & pred_latlong$long_dist <0, 1, ifelse(abs(pred_latlong$long_dist)>180 & pred_latlong$long_dist >0, 3, 2))
pred_latlong$distance <- distHaversine(p1 = pred_latlong[, c("source_longitude", "source_latitude")], p2=pred_latlong[, c("island_longitude", "island_latitude")])

#Average distance (in km) between source and island
mean(pred_latlong$distance)/1000

p2 <- ggplot()+
  geom_sf(data=world) +
  xlab("Longitude") + ylab("Latitude") +
  geom_point(data=pred_latlong, aes(y=island_latitude, x=island_longitude), color="blue")+
  geom_point(data=pred_latlong, aes(y=source_latitude, x=source_longitude))+
  theme_cowplot()+
  theme(legend.position="none")+
  geom_segment(data=subset(pred_latlong, direction==2), aes(y=source_latitude, yend=island_latitude, x=source_longitude, xend=island_longitude), linewidth=0.1)+
 geom_segment(data=subset(pred_latlong, direction==1), aes(y=source_latitude, yend=island_latitude, x=source_longitude, xend=-Inf), linewidth=0.1)+
  geom_segment(data=subset(pred_latlong, direction==1), aes(y=source_latitude, yend=island_latitude, x=Inf, xend=island_longitude), linewidth=0.1)+
 geom_segment(data=subset(pred_latlong, direction==3), aes(y=source_latitude, yend=island_latitude, x=source_longitude, xend=Inf), linewidth=0.1)+
 geom_segment(data=subset(pred_latlong, direction==3), aes(y=source_latitude, yend=island_latitude, x=-Inf, xend=island_longitude), linewidth=0.1)
p2
save_plot("top_latlong_map_ants.pdf", p2)

save(pred_latlong, file = "data/konig/gdm_pred_latlong.RData")

#############################
### Aggregate predictions ###
#############ä###############

# Since the predictions are for grid cells, not the actual mainland units, we have to aggregate the predicted values
# Get ant polygons per mainland entity
mainland_gridcells_ants = lapply(Bentity2$BENTITY2_N, function(ml_ID){
  print(ml_ID)
  tmp_entity = Bentity2[which(Bentity2$BENTITY2_N == ml_ID),]
  tryCatch({
    tmp_grid = raster::intersect(grid, tmp_entity)
    list(ID = tmp_grid@data$ID,
         area = raster::area(tmp_grid)/1000000)},
    error = function(e){NA})
})

names(mainland_gridcells_ants) = Bentity2$BENTITY2_N

# Prepare source matrix
ant_source_matrix = matrix(NA, nrow = length(islands), ncol = length(Bentity2$BENTITY2_N), byrow = T)
rownames(ant_source_matrix) = islands
colnames(ant_source_matrix) = Bentity2$BENTITY2_N

# Fill source matrix
for(island in islands){
  print(island)
  island_richness = length(which(gabi_all$entity_ID == island))
  tmp_sources = sapply(as.character(Bentity2$BENTITY2_N), function(x){
    tmp_grid = mainland_gridcells_ants[[x]]
    tmp_grid[[1]] <- as.numeric(tmp_grid[[1]] + 100000)
    if(is.na(any(tmp_grid[[1]])) || length(tmp_grid[[1]]) == 0){return(NA)}
    grid_subset = which(tmp_grid[[1]] %in% colnames(gdm_pred_mat))
    grid_IDs = tmp_grid[[1]][grid_subset]
    grid_areas = tmp_grid[[2]][grid_subset]
    if(is.na(any(grid_IDs)) || length(grid_IDs) == 0){return(NA)}
    grid_index = which(colnames(gdm_pred_mat) %in% grid_IDs)
    grid_index = grid_index[which(!is.na(grid_index))]
    grid_values = gdm_pred_mat[paste(island), grid_index]
    p = 1 - weighted.mean(grid_values, grid_areas)
  })
  names(tmp_sources) = Bentity2$BENTITY2_N
  ant_source_matrix[paste(island),] = tmp_sources # write in matrix
}

save(ant_source_matrix, file = "data/konig/ant_source_matrix.RData")

##############
### Output ###
##############

# Extract top 5 highest predicted mainlands from GDM per island
pairs_pred <- data.frame(island=NULL, mainland_pred=NULL, rank=NULL)

ant_source_matrix <- ant_source_matrix[, colnames(ant_source_matrix) %in% meta_mainlands$geo_entity]

for(i in 1:length(rownames(ant_source_matrix))){
  tmp <- as.data.frame(ant_source_matrix[i, ])
  colnames(tmp) <- "predict"
  topfive <- tmp %>% arrange(-predict) %>% head(5)
  topfive$island <- rownames(ant_source_matrix)[[i]]
  topfive$mainland_pred <- rownames(topfive)
  topfive$rank <- rank(-topfive$predict)
  pairs_pred <- rbind(pairs_pred, topfive)
}

#add in island metadata
pairs_pred <- merge(pairs_pred, meta_islands, by.x="island", by.y="entity_ID", all.x=T, all.y=F)
#add in mainland metadata
pairs_pred <- merge(pairs_pred, meta_mainlands, by.x="mainland_pred", by.y="geo_entity", all.x=T, all.y=F)
pairs_pred$island <- as.numeric(pairs_pred$island)

save(pairs_pred, file = "data/konig/pairs_pred_ants.RData")

############################################################
### Visualize predicted source mainlands for each island ###
############################################################

#Make map of most likely source mainland for each island
pairs_pred$rank <- as.factor(pairs_pred$rank)

#A little code to prevent lines from crossing between the South Pacific on the left and Australia on the right
pairs_pred$long_dist <- pairs_pred$longitude.x-pairs_pred$longitude.y
pairs_pred$direction <- ifelse(abs(pairs_pred$long_dist)>180 & pairs_pred$long_dist <0, 1, ifelse(abs(pairs_pred$long_dist)>180 & pairs_pred$long_dist >0, 3, 2))

#Calculate distances betweens sources and islands
pairs_pred$distance <- distHaversine(p1 = pairs_pred[, c("longitude.x", "latitude.x")], p2=pairs_pred[, c("longitude.y", "latitude.y")])

#Average distance (in km) between source and island
mean(pairs_pred$distance)/1000

p3 <- ggplot()+
  geom_sf(data=world) +
  xlab("Longitude") + ylab("Latitude") +
  geom_point(data=pairs_pred, aes(y=latitude.x, x=longitude.x), color="blue")+
  geom_point(data=pairs_pred, aes(y=latitude.y, x=longitude.y, color=rank))+
  theme_cowplot()+
  labs(color="Rank")+
  geom_segment(data=subset(pairs_pred, direction==2), aes(y=latitude.x, yend=latitude.y, x=longitude.x, xend=longitude.y), linewidth=0.1)+
  geom_segment(data=subset(pairs_pred, direction==1), aes(y=latitude.x, yend=latitude.y, x=longitude.x, xend=-Inf), linewidth=0.1)+
  geom_segment(data=subset(pairs_pred, direction==1), aes(y=latitude.x, yend=latitude.y, x=Inf, xend=longitude.y), linewidth=0.1)+
  geom_segment(data=subset(pairs_pred, direction==3), aes(y=latitude.x, yend=latitude.y, x=longitude.x, xend=Inf), linewidth=0.1)+
  geom_segment(data=subset(pairs_pred, direction==3), aes(y=latitude.x, yend=latitude.y, x=-Inf, xend=longitude.y), linewidth=0.1)
p3
save_plot("topfive_map_ants.pdf", p3)

p4 <- ggplot()+
  geom_sf(data=world) +
  xlab("Longitude") + ylab("Latitude") +
  geom_point(data=subset(pairs_pred, rank == 1), aes(y=latitude.x, x=longitude.x), color="blue")+
  geom_point(data=subset(pairs_pred, rank == 1), aes(y=latitude.y, x=longitude.y, color=rank))+
  theme_cowplot()+
  theme(legend.position="none")+
  geom_segment(data=subset(pairs_pred, direction==2 & rank == 1), aes(y=latitude.x, yend=latitude.y, x=longitude.x, xend=longitude.y), linewidth=0.1)+
  geom_segment(data=subset(pairs_pred, direction==1 & rank == 1), aes(y=latitude.x, yend=latitude.y, x=longitude.x, xend=-Inf), linewidth=0.1)+
  geom_segment(data=subset(pairs_pred, direction==1 & rank == 1), aes(y=latitude.x, yend=latitude.y, x=Inf, xend=longitude.y), linewidth=0.1)+
  geom_segment(data=subset(pairs_pred, direction==3 & rank == 1), aes(y=latitude.x, yend=latitude.y, x=longitude.x, xend=Inf), linewidth=0.1)+
  geom_segment(data=subset(pairs_pred, direction==3 & rank == 1), aes(y=latitude.x, yend=latitude.y, x=-Inf, xend=longitude.y), linewidth=0.1)
p4
save_plot("top_map_ants.pdf", p4)
