################################################################################
#### Project: Festuca metabolites
#### Title:   Function | Load climate niche data
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    21 February 2022
#### ---------------------------------------------------------------------------

load_climate <- function(rawPlants){
  ## Load raw data ----
  # data tables
  dataFrb <- fread("./data/climrubra.csv", data.table = F)
  dataDEM <- fread("./data/Taxon_Predictors_Ranges.csv", data.table = F)
  # raster files
  tifFiles <- list.files(
    "./data/rasters",
    pattern = ".tif",
    full.names = T
  )
  tifData <- lapply(tifFiles, raster)
  names(tifData) <- substr(tifFiles, 16, nchar(tifFiles) - 4)
  ## Process raster data ----
  # reproject coordinates onto Swiss grid
  for(i in 1:length(tifData)){
    tifData[[i]] <- projectRaster(
      tifData[[i]], 
      crs = "+init=epsg:21781"
    )
    cat("i\n")
  }
  # stack data together
  tifStack <- stack(tifData)
  # get xy coordinates
  xy <- select(dataFrb, Longitude, Latitude)
  coordinates(xy) <- ~ Latitude + Longitude
  # get data
  climFrb <- raster::extract(tifStack, xy)
  ## Format DEM data ----
  # split names
  splitNames <- str_split(dataDEM$taxon, "_") %>%
    lapply(., function(x){x[1:2]}) %>%
    do.call(rbind, .)
  # format full dataset
  dataDEMclean <- dataDEM %>%
    # create species identifiers
    mutate(genus = splitNames[, 1]) %>%
    mutate(species = splitNames[, 2]) %>%
    mutate(genusSpecies = paste(genus, species)) %>%
    # select columns of interest
    select(genusSpecies, genus, species, DEM_Q05:SoilR_Q95) %>%
    # filter for festuca only
    filter(genusSpecies %in% rawPlants$cross$species) %>%
    # summarise data at species level
    group_by(genusSpecies, genus, species) %>%
    summarise(across(everything(), mean)) %>%
    # filter for columns of interest
    select(-contains("Q05")) %>%
    select(-contains("Q95")) %>%
    select(-contains("_Median")) %>%
    # rename remaining columns
    select(
      genusSpecies:species,
      dem = DEM_Mean,
      ndvi = ndviMEAN_Mean,
      srad = sradyy_Mean,
      tave = taveyy_Mean,
      prec = precyy_Mean,
      forest = ForestQ25_Mean,
      TPI = TPI_Mean,
      TRI = TRI_Mean,
      TWI = TWI_Mean,
      soilF = SoilF_Mean,
      soilR = SoilR_Mean
    ) %>%
    ungroup %>%
    as.data.frame
  ## Format F. rubra data ----
  dataFrbClean <- dataFrb %>%
    # change altitude to lower case
    mutate(elevDisc = tolower(Alt)) %>%
    bind_cols(., climFrb) %>%
    # select columns
    select(
      transect = Transect,
      # note switched lon/lat!
      lat = Longitude,
      lon = Latitude,
      elevCont = Altitude,
      elevDisc,
      ndvi = ClimHerb_ndviMEAN,
      srad = Daymet_sradyy,
      tave = Daymet_taveyy,
      prec = Daymet_precyy,
      forest = MoGLI_ForestQ25,
      TPI = R_TPI,
      TRI = R_TRI,
      TWI = SAGA_TWI,
      soilF = SPEEDMIND_SoilF,
      soilR = SPEEDMIND_SoilR
    )
  ## Collate and return ----
  climOut <- list(
    rubra = dataFrbClean,
    cross = dataDEMclean
  )
  return(climOut)
}
