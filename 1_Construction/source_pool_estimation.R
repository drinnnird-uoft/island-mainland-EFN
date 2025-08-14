#Based on Konig et al. (2020) in Ecography
#See https://datadryad.org/dataset/doi:10.5061/dryad.6wwpzgmwm


library(tidyverse)
library(vegan)
library(gdm)
library(geosphere)
library(raster)
library(cowplot)
library("rnaturalearth")
library("rnaturalearthdata")
library(ggspatial)
library(segmented)

#Load checklist data
load("data/konig/checklists.RData")
#Load meta-data for 178 islands and 735 mainlands
load("data/konig/geoentities_meta.RData")
#Load entity IDs that have floral checklists for main analysis, to filter set of possible source mainlands
load("GIFT_UniqueEntities.RData") # 1771 unique geo entities

#Filter geo entities by unique entities with floral checklists
geoentities_meta <- geoentities_meta %>% 
  filter(entity_ID %in% unique_geo_entities)

#Check overlap between checkliss and islands and mainlands with metadata
length(checklists$entity_ID) #1669603 rows
length(geoentities_meta) #713 entities
sum(checklists$entity_ID %in% geoentities_meta$entity_ID) #Most checklist entity IDs are in the metadata; 1222440
sum(geoentities_meta$entity_ID %in% checklists$entity_ID) #All 713 of the metadata IDs are in the checklist


#I think the work_ID is a unique identifier for each species


###########################
### Calculate distances ###
###########################
# 1. Compositional dissimilarity (Species level)
spec_tab = unclass(table(checklists[,c("entity_ID", "work_ID")])) # look at species composition
#spec_tab is a matrix of sites (1036 entities in the checklists dataset) by species (240961 work_IDs)
spec_dist = betadiver(spec_tab, "sim") #   Turnover / Simpsons Beta diversity
spec_dist = as.matrix(spec_dist) # Convert to matrix


# 2. Geographical distance
geo_dist = round(distm(as.matrix(geoentities_meta[, c("longitude", "latitude")]), fun = distHaversine)/1000,1)+0.1
colnames(geo_dist) = geoentities_meta$entity_ID
rownames(geo_dist) = geoentities_meta$entity_ID


save(list = c("geo_dist", "spec_dist"), file = "data/konig/distances.RData")
gc() #Garbage collection


############################################
###  Generalized dissimilarity modeling  ###
############################################
# We fit a model to reflect the default behaviour of floristic similarity on the mainland 
# using distance + climatic variables. Given new data, we can predict the expected number of shared species
# relative to a global equal area grid in order to localize potential source pools for islands. 
spec_dist = cbind("entity_ID" = as.numeric(rownames(spec_dist)), spec_dist)
geo_dist = cbind("entity_ID" = as.numeric(rownames(geo_dist)), geo_dist)


#This just makes a list of mainland entity IDs and a list of island entity ID
mainlands = geoentities_meta$entity_ID[geoentities_meta$entity_class %in% c("Mainland", "Island/Mainland")]
islands = geoentities_meta$entity_ID[geoentities_meta$entity_class == "Island"]


### 1. Model fitting
spec_dist_ml = spec_dist[paste(mainlands), c("entity_ID", paste(mainlands))]
geo_dist_ml = geo_dist[paste(mainlands), c("entity_ID", paste(mainlands))]
env_data_ml = geoentities_meta %>% filter(entity_ID %in% mainlands) %>% dplyr::select(entity_ID, longitude, latitude, T_mean, T_var, P_mean, P_var)


#This fits a model using mainland data
gdm_tab_ml = formatsitepair(bioData = spec_dist_ml, bioFormat = 3, siteColumn = "entity_ID", abundance = F, XColumn = "longitude", YColumn = "latitude", predData = env_data_ml, distPreds = list("Geo_Distance" = geo_dist_ml))
gdm_fit_ml = gdm(gdm_tab_ml, geo = F)
summary(gdm_fit_ml)


#Save the model
save(gdm_fit_ml, file = "data/konig/gdm_fit_ml.RData")


### 2. Prepare prediction Data
# load global equal area grid and prepared environmental variables for each grid cell
load(file = "data/konig/grid.RData")
load(file = "data/konig/env_grid.RData")


# Combine environmental data of grid cells and islands 
env_data_is = geoentities_meta %>% filter(entity_ID %in% islands) %>% dplyr::select(entity_ID, longitude, latitude, T_mean, T_var, P_mean, P_var)
#Make entity ID numeric, as it no longer can be a string
env_grid$entity_ID <- gsub("cell_", "", env_grid$entity_ID)


# To distinguish 'real' entities from grid cell entities, add 100000 to grid cell entity IDs
#Any value over 100000 is from the grid
env_grid$entity_ID <- as.numeric(env_grid$entity_ID)+100000
env_grid = rbind(env_data_is, env_grid)


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
save(gdm_tab_pred, file = "data/konig/gdm_tab_pred.RData")


### 3. Predict
#Changed from predict.gdm to predict, but seems to work fine
gdm_pred = predict(gdm_fit_ml, gdm_tab_pred, geo = F, time = F)
gdm_pred_mat = matrix(NA, nrow = nrow(env_grid), ncol = nrow(env_grid))
gdm_pred_mat[which(lower.tri(gdm_pred_mat))] = gdm_pred
gdm_pred_mat = as.matrix(as.dist(gdm_pred_mat))


entity_names = sort(env_grid$entity_ID)
colnames(gdm_pred_mat) = entity_names
rownames(gdm_pred_mat) = entity_names

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

save(gdm_pred_mat, file = "data/konig/gdm_pred_mat.RData")
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

#World map base layer
world <- ne_countries(scale = "medium", returnclass = "sf")

p1 <- ggplot()+
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
p1
save_plot("top_latlong_map.pdf", p1)

save(pred_latlong, file = "data/konig/gdm_pred_latlong.RData")

#############################
### Aggregate predictions ###
#############ä###############
load("data/konig/geoentities.RData")

# Since the predictions are for grid cells, not the actual mainland units, we have to aggregate the predicted values
# Get cell IDs per mainland entity
mainland_gridcells = lapply(mainlands, function(ml_ID){
  print(ml_ID)
  tmp_entity = geoentities[which(geoentities@data$entity_ID == ml_ID),]
  tryCatch({
    tmp_grid = raster::intersect(grid, tmp_entity)
    list(ID = tmp_grid@data$ID,
         area = raster::area(tmp_grid)/1000000)},
    error = function(e){NA})
})


names(mainland_gridcells) = mainlands


# Prepare source matrix
source_matrix = matrix(NA, nrow = length(islands), ncol = length(mainlands), byrow = T)
rownames(source_matrix) = islands
colnames(source_matrix) = mainlands


# Fill source matrix
for(island in islands){
  print(island)
  island_richness = length(which(checklists$entity_ID == island))
  tmp_sources = sapply(mainlands, function(x){
    tmp_grid = mainland_gridcells[[paste(x)]]
    tmp_grid$ID <- as.numeric(tmp_grid$ID + 100000)
    if(is.na(any(tmp_grid$ID)) || length(tmp_grid$ID) == 0){return(NA)}
    grid_subset = which(tmp_grid$ID %in% colnames(gdm_pred_mat))
    grid_IDs = tmp_grid$ID[grid_subset]
    grid_areas = tmp_grid$area[grid_subset]
    if(is.na(any(grid_IDs)) || length(grid_IDs) == 0){return(NA)}
    grid_index = which(colnames(gdm_pred_mat) %in% grid_IDs)
    grid_index = grid_index[which(!is.na(grid_index))]
    grid_values = gdm_pred_mat[paste(island), grid_index]
    p = 1 - weighted.mean(grid_values, grid_areas)
  })
  names(tmp_sources) = mainlands
  
  # Remove all mainland units beyond "knee" (see Supporting Text 2 for justification)
  y = 1- sort(tmp_sources, decreasing = T)
  y_norm = (y-min(y)) / (max(y)-min(y)) # Normalize y to (0,1) to make slopes comparable
  spline_fit = smooth.spline(1:length(y), y_norm, spar = 0.75) # Best for this purpose, because smooth and cont. increasing
  spline_prime = diff(spline_fit$y) # Derivative, i.e. what is the slope of 'spline_fit'?
  spline_runlength = rle(spline_prime < 0.001) # Consecutive slope values smaller/greater than 0.001
  if(tail(spline_runlength$values, n = 1) == T){ # Is the last run < 0.001, i.e. the function flattens out?
    cutoff = length(y) - tail(spline_runlength$lengths, n = 1) # cutoff value
  } else {
    cutoff = length(y) - sum(tail(spline_runlength$lengths, n = 2))
  }
  tmp_sources[names(y[cutoff:length(y)])] = 0 # set all mainlands beyond cutoff to zero
  source_matrix[paste(island),] = tmp_sources # write in matrix
}


save(source_matrix, file = "data/konig/source_matrix.RData")
#load("data/konig/source_matrix.RData")


#######################
### Validate output ###
#######################


#One thing I want to check is whether the GDM models make a different prediction about source mainlands
#than simply examining species turnover between islands and mainlands (e.g., low turnove=likely source mainland)
#This info is contained in the spec_dist matrix

# Extract top 5 mainlands with most similar species composition as each island from spec_dist matrix
pairs_beta <- data.frame(island=NULL, mainland_beta=NULL)

for(i in 1:dim(spec_dist)[[1]]){
  tmp <- as.data.frame(spec_dist[i, ])
  colnames(tmp) <- "beta"
  topfive <- tmp %>% arrange(beta) %>% head(6)
  topfive$island <- tmp[rownames(tmp)=="entity_ID", ]
  topfive$mainland_beta <- rownames(topfive)
  pairs_beta <- rbind(pairs_beta, topfive[-1, ])
}

# Extract top 5 highest predicted mainlands from GDM per island
pairs_pred <- data.frame(island=NULL, mainland_pred=NULL, rank=NULL)

for(i in 1:length(rownames(source_matrix))){
  tmp <- as.data.frame(source_matrix[i, ])
  colnames(tmp) <- "predict"
  topfive <- tmp %>% arrange(-predict) %>% head(5)
  topfive$island <- rownames(source_matrix)[[i]]
  topfive$mainland_pred <- rownames(topfive)
  topfive$rank <- rank(-topfive$predict)
  pairs_pred <- rbind(pairs_pred, topfive)
}

all_pairs <- merge(pairs_pred, pairs_beta, by="island", all.x=T, all.y=F)
all_pairs$match <- all_pairs$mainland_pred == all_pairs$mainland_beta


#What percentage of source mainlands are the same using just beta-diversity compared to the GDM model
sum(all_pairs$match)/178*100 #Twenty-one percent


#Okay, keep going with the GDM predicted mainlands
#add in island metadata
pairs_pred <- merge(pairs_pred, geoentities_meta, by.x="island", by.y="entity_ID", all.x=T, all.y=F)
#add in mainland metadata
pairs_pred <- merge(pairs_pred, geoentities_meta, by.x="mainland_pred", by.y="entity_ID", all.x=T, all.y=F)
pairs_pred$island <- as.numeric(pairs_pred$island)

saveRDS(pairs_pred, "data/konig/final_predicted_mainland_island_pairs.Rds")


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

p2 <- ggplot()+
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
p2
save_plot("topfive_map.pdf", p2)

p3 <- ggplot()+
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
p3
save_plot("top_map.pdf", p3)

#How many mainlands should we use for each island? 
multi_plot <- ggplot()
for(island in islands){
  tmp <- as.data.frame(source_matrix[which(rownames(source_matrix) == island), ])
  tmp$mainland <- rownames(tmp)
  colnames(tmp) <- c("island", "mainland")
  tmp$mainland_rank <- length(tmp$mainland)-rank(tmp$island)
  multi_plot <- multi_plot+geom_line(data=tmp, aes(y=island, x=mainland_rank, group=1), linewidth=0.2, alpha=0.5)+theme_cowplot()+xlab("Mainland rank")+ylab("Probability of being source mainland")
}
multi_plot <- multi_plot+scale_x_continuous(limits=c(0,25))
save_plot("decay_with_rank.pdf", multi_plot)
#Maybe 5 is a good number

