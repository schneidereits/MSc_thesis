#### validation
library(tidyverse)
library(terra)
library(tidyterra)
library(sf)

# import wetland mask for setting projections
wetland_mask <- rast("/data/Dagobah/greengrass/schnesha/thesis/wetland_masks/mosaic/wetland_mask.vrt")

# load meta data of VHI data
BB_metadata <- vect("/data/Dagobah/greengrass/schnesha/thesis/05_validation/BB_2020-2022_DOP_meta_data/creationdate.shp") |> 
  rename(date_validation=creationda) |> 
  select(date_validation) |> 
  mutate(date_validation = as.Date(date_validation)) |>  project(wetland_mask)

plot(BB_metadata)


# add sample point data
sample_points_wide <- readRDS( file = "/data/Dagobah/greengrass/schnesha/thesis/feature_space/sample_points_df_regular_wide.rds") |> 
  filter(date %in% unique(BB_metadata$date_validation)) |> slice_sample(n=1000000)

# convert to vect
sample_points_vect <- vect(sample_points_wide, geom=c("x","y"))
crs(sample_points_vect)  <- "epsg:3035"
sample_points_vect <- sample_points_vect |> project(wetland_mask) |> mask(BB_metadata)

plot(sample_points_vect,add=T)

# intersect sample points and meta data
# CRITICAL ERROR THAT OCCURES FOR NO APPERNT REASON.
# foo <- terra::intersect(BB_metadata, sample_points_vect)

sample_points <- st_intersection(sf::st_as_sf(BB_metadata), sf::st_as_sf(sample_points_vect)) |>
  vect() |> 
  filter(date >= (date_validation - 5) & date <= (date_validation + 5)) |> 
  as.data.frame(geom="XY")

  
  # set bin with for feature space stratificaiton
  NDVI_range = c(seq(-1,1,0.25))
  SWIR_range =c(seq(0,1,0.25))
  bins <- expand_grid(NDVI_range, SWIR_range)
  
  
  validation <- c()
  for (i in 1:(length(NDVI_range)-1)) {
    for (j in 1:(length(SWIR_range)-1)) {
      
      temp <- sample_points |> 
        filter(NDVI>NDVI_range[i],
               NDVI<NDVI_range[i+1],
               SWIR_ratio>SWIR_range[j],
               SWIR_ratio<SWIR_range[j+1]) 
      
      # limit water points in validation  
      if (nrow(temp)>10 & NDVI_range[i+1]<0) {
        temp <- temp |>
          slice_sample(n=6)
      }
      
      # reduce size if to many points are in a feature space bin
      if (nrow(temp)>20) {
        temp <- temp |>
          slice_sample(n=20)
      }
      
      
      # bind to final output df
      validation <- plyr::rbind.fill(validation,temp)
      
    }
  }
  
  validation |> 
    group_by(month_year, Tile_ID) |> 
    count()
  
  validation |> 
    group_by(year) |> 
    count()
  
  
  ggplot(validation, aes(NDVI, SWIR_ratio)) +
    geom_point() +
    scale_fill_continuous(type = "viridis") +
    theme_minimal() 
  
  
  ### save to vector 
  
  
  xy <- validation |> select(x,y) 
  validation_points <- validation |> select(x, y, ID, Tile_ID, date, month, year, month_year, doy, UID ) |> 
    mutate(NPV=-999,
           PV=-999,
           soil=-999,
           water=-999)
  
  library(terra)
  library(sp)
  validation_points_vect <- SpatialPointsDataFrame(coords = xy, 
                                                   data = subset(validation_points),
                                                   proj4string = CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs +type=crs"))  |> 
    vect()
  
  plot(validation_points_vect)
  
  writeVector(validation_points_vect, 
              "/data/Dagobah/greengrass/schnesha/thesis/05_validation/validation_points.gpkg",
              filetype="GPKG",
              options = NULL,
              overwrite=TRUE)
  
### reimport validated points 
  
validated_points_vect <- vect("/data/Dagobah/greengrass/schnesha/thesis/05_validation/validation_points_validated.gpkg")
validated_points <-
 # as.data.frame(validated_points_vect) |> 
  validated_points_vect |> 
    filter(!(NPV == -999 & PV == -999 & soil == -999 & water == -999)) |> 
    mutate(across(c(NPV, PV, soil, water), ~ifelse(. == -999, 0, .)),
           total_cover = NPV + PV + soil + water,
           date= as.numeric(date),
           date = as.Date(as.POSIXct(date, origin="1970-01-01")))
  

FC_parameters <- as.data.frame(validated_points) %>% distinct(Tile_ID, date)

FC_files <- c()
validation_predition_df <- c()
for (i in 1:nrow(FC_parameters)) {
  
  temp_validation <- validated_points |> 
    filter(Tile_ID==FC_parameters[i,1],
           date== (FC_parameters[i,2]))
  
  
  temp_prediction <- list.files(paste0("/data/Dagobah/greengrass/schnesha/thesis/04_fraction_mapping_ts/05_predictions/",
                           FC_parameters[i,1]), full.names = T)
  
  date <- str_replace_all(FC_parameters[i,2],"-","")
  temp_prediction <- temp_prediction[grepl(date, temp_prediction)]
  
  FC_files[i] <- temp_prediction
  
  validation_predition_df <- plyr::rbind.fill(validation_predition_df,
                                         as.data.frame(terra::extract(rast(temp_prediction), temp_validation,  bind=T)))

}

validation_df <- validation_predition_df |> 
  rename(NPV_prediction = MODEL_CLASS_001_ITERATION_001,
         PV_prediction = MODEL_CLASS_002_ITERATION_001,
         soil_prediction = MODEL_CLASS_003_ITERATION_001,
         water_prediction = MODEL_CLASS_004_ITERATION_001) |> 
  mutate(across(c(NPV_prediction:water_prediction), ~./100))


### PV
(p1 <- ggplot(validation_df, aes(x=PV_prediction, y=PV)) +
  geom_point()+
  xlim(0,100) +
  ylim(0,100) +
  scale_color_viridis()+
  geom_abline(aes(slope=1, intercept=0))+
  geom_smooth(method="lm", color="olivedrab")+
  annotate("text",
           x=1,
           y=100,
           label= paste0("R2 =",round((cor(validation_df$PV, validation_df$PV_prediction, use="complete.obs")^2),2)),
           size=4.5,
           hjust=0)+
  annotate("text",
           x= 1, y=95,
           hjust=0, size=4.5,
           label= paste0("RMSE ==",round(sqrt(mean((validation_df$PV- validation_df$PV_prediction)^2, na.rm=TRUE)),2)),
           parse=TRUE)+
  annotate("text",
           x=1, y= 90, 
           size=4.5, hjust=0,
           label=paste0("Bias = ", round((mean(validation_df$PV, na.rm=TRUE) - mean(validation_df$PV_prediction,na.rm=TRUE)),2)))+
  annotate("text",
           x=1, y=85,
           size=4.5,hjust=0,
           label=paste0("MAE = ", round(mean(abs(validation_df$PV - validation_df$PV_prediction)),2))) +
  theme_classic()
)

# NPV
(p2 <- ggplot(validation_df, aes(x=NPV_prediction, y=NPV)) +
  geom_point()+
  xlim(0,100) +
  ylim(0,100) +
  scale_color_viridis()+
  geom_abline(aes(slope=1, intercept=0))+
  geom_smooth(method="lm", color="peru")+
  annotate("text",
           x=1,
           y=100,
           label= paste0("R2 =",round((cor(validation_df$NPV, validation_df$NPV_prediction, use="complete.obs")^2),2)),
           size=4.5,
           hjust=0)+
  annotate("text",
           x= 1, y=95,
           hjust=0, size=4.5,
           label= paste0("RMSE ==",round(sqrt(mean((validation_df$NPV- validation_df$NPV_prediction)^2, na.rm=TRUE)),2)),
           parse=TRUE)+
  annotate("text",
           x=1, y= 90, 
           size=4.5, hjust=0,
           label=paste0("Bias = ", round((mean(validation_df$NPV, na.rm=TRUE) - mean(validation_df$NPV_prediction,na.rm=TRUE)),2)))+
  annotate("text",
           x=1, y=85,
           size=4.5,hjust=0,
           label=paste0("MAE = ", round(mean(abs(validation_df$NPV - validation_df$NPV_prediction)),2))) +
  theme_classic()
)

# Soil
(p3 <-ggplot(validation_df, aes(x=soil_prediction, y=soil)) +
  geom_point()+
  xlim(0,100) +
  ylim(0,100) +
  scale_color_viridis()+
  geom_abline(aes(slope=1, intercept=0))+
  geom_smooth(method="lm", color="darkred")+
  annotate("text",
           x=1,
           y=100,
           label= paste0("R2 =",round((cor(validation_df$soil, validation_df$soil_prediction, use="complete.obs")^2),2)),
           size=4.5,
           hjust=0)+
  annotate("text",
           x= 1, y=95,
           hjust=0, size=4.5,
           label= paste0("RMSE ==",round(sqrt(mean((validation_df$soil- validation_df$soil_prediction)^2, na.rm=TRUE)),2)),
           parse=TRUE)+
  annotate("text",
           x=1, y= 90, 
           size=4.5, hjust=0,
           label=paste0("Bias = ", round((mean(validation_df$soil, na.rm=TRUE) - mean(validation_df$soil_prediction,na.rm=TRUE)),2)))+
  annotate("text",
           x=1, y=85,
           size=4.5,hjust=0,
           label=paste0("MAE = ", round(mean(abs(validation_df$soil - validation_df$soil_prediction)),2))) +
  theme_classic()
)
# Water
(p4 <- ggplot(validation_df, aes(x=water_prediction, y=water)) +
  geom_point()+
  xlim(0,100) +
  ylim(0,100) +
  scale_color_viridis()+
  geom_abline(aes(slope=1, intercept=0))+
  geom_smooth(method="lm", color="steelblue")+
  annotate("text",
           x=1,
           y=100,
           label= paste0("R2 =",round((cor(validation_df$water, validation_df$water_prediction, use="complete.obs")^2),2)),
           size=4.5,
           hjust=0)+
  annotate("text",
           x= 1, y=95,
           hjust=0, size=4.5,
           label= paste0("RMSE ==",round(sqrt(mean((validation_df$water- validation_df$water_prediction)^2, na.rm=TRUE)),2)),
           parse=TRUE)+
  annotate("text",
           x=1, y= 90, 
           size=4.5, hjust=0,
           label=paste0("Bias = ", round((mean(validation_df$water, na.rm=TRUE) - mean(validation_df$water_prediction,na.rm=TRUE)),2)))+
  annotate("text",
           x=1, y=85,
           size=4.5,hjust=0,
           label=paste0("MAE = ", round(mean(abs(validation_df$water - validation_df$water_prediction)),2))) +
  theme_classic()
)

(p1 +p2 +p3 +p4)


# old sample point generation method
sample_points_wide <- readRDS( file = "/data/Dagobah/greengrass/schnesha/thesis/feature_space/sample_points_df_regular_wide.rds") %>% slice_sample(n=100000)


ggplot(sample_points_wide, aes(NDVI, SWIR_ratio)) +
  stat_bin_hex(bins = 200) +
  scale_fill_continuous(type = "viridis") +
  theme_minimal() 

# set bin with for feature space stratificaiton
NDVI_range = c(seq(-1,1,0.25))
SWIR_range =c(seq(0,1,0.25))
bins <- expand_grid(NDVI_range, SWIR_range)


validation <- c()
for (i in 1:(length(NDVI_range)-1)) {
  for (j in 1:(length(SWIR_range)-1)) {
    
    temp <- sample_points_wide |> 
      filter(NDVI>NDVI_range[i],
             NDVI<NDVI_range[i+1],
             SWIR_ratio>SWIR_range[j],
             SWIR_ratio<SWIR_range[j+1]) |> 
      # change to doy (n-7):(n+7)
      filter(Tile_ID== "X0070_Y0038" & year=="2018" & doy %in% c((156-5):(156+5)) |
               Tile_ID== "X0068_Y0037" & year=="2018" & doy %in% c((156-5):(156+5)) |
               Tile_ID== "X0068_Y0037" & year=="2020" & doy %in% c((225-5):(225+5)) |
               Tile_ID== "X0068_Y0037" & year=="2022" & doy %in% c((101-5):(101+5)) |
               Tile_ID== "X0068_Y0038" & year=="2018" & doy %in% c((156-5):(156+5)) |
               Tile_ID== "X0068_Y0038" & year=="2020" & doy %in% c((225-5):(225+5)) |
               Tile_ID== "X0068_Y0038" & year=="2022" & doy %in% c((101-5):(101+5)) |
               Tile_ID== "X0068_Y0042" & year=="2019" & doy %in% c((110-5):(110+5)) |
               Tile_ID== "X0068_Y0042" & year=="2020" & doy %in% c((109-5):(109+5)) |
               Tile_ID== "X0068_Y0042" & year=="2021" & doy %in% c((151-5):(151+5)) |
               Tile_ID== "X0071_Y0045" & year=="2017" & doy %in% c((152-5):(152+5)) |
               Tile_ID== "X0071_Y0045" & year=="2020" & doy %in% c((264-5):(264+5)) |
               Tile_ID== "X0071_Y0045" & year=="2021" & doy %in% c((251-5):(251+5)))
    
    # limit water points in validation  
    if (nrow(temp)>10 & NDVI_range[i+1]<0) {
      temp <- temp |>
        slice_sample(n=3)
    }
    
    # reduce size if to many points are in a feature space bin
    if (nrow(temp)>10) {
      temp <- temp |>
        slice_sample(n=10)
    }
    
    
    # bind to final output df
    validation <- plyr::rbind.fill(validation,temp)
    
  }
}

validation |> 
  group_by(month_year, Tile_ID) |> 
  count()

validation |> 
  group_by(year) |> 
  count()


ggplot(validation, aes(NDVI, SWIR_ratio)) +
  geom_point() +
  scale_fill_continuous(type = "viridis") +
  theme_minimal() 


### save to vector 


xy <- validation |> select(x,y) 
validation_points <- validation |> select(x, y, ID, Tile_ID, date, month, year, month_year, doy, UID ) |> 
  mutate(NPV=-999,
         PV=-999,
         soil=-999,
         water=-999)

library(terra)
library(sp)
validation_points_vect <- SpatialPointsDataFrame(coords = xy, 
                                                 data = subset(validation_points),
                                                 proj4string = CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs +type=crs"))  |> 
  vect()

plot(validation_points_vect)

writeVector(validation_points_vect, 
            "/data/Dagobah/greengrass/schnesha/thesis/05_validation/validation_points.gpkg",
            filetype="GPKG",
            options = NULL,
            overwrite=TRUE)



