################################################################################
#### Project: Festuca metabolites
#### Title:   Function | Load climate niche data
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    21 February 2022
#### ---------------------------------------------------------------------------

load_climate <- function(rawPlants){
  ## Load raw data ----
  dataDEM <- fread("./data/Taxon_Predictors_Ranges.csv", data.table = F)
  dataFrb <- fread("./data/climrubra.csv", data.table = F)
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
    summarise(across(everything(), median)) %>%
    ungroup %>%
    as.data.frame
  ## Format F. rubra data ----
  dataFrbClean <- dataFrb %>%
    # change altitude to lower case
    mutate(elevDisc = tolower(Alt)) %>%
    # select columns
    select(
      transect = Transect,
      elevCont = Altitude,
      elevDisc,
      lon = Longitude,
      lat = Latitude,
      vec_srad59:vec_sfroyy
    ) %>%
    # rename climate data to remove precursory dataset name
    rename_with(function(x){substr(x, 5, nchar(x))}, starts_with("vec_"))
  ## Collate and return ----
  climOut <- list(
    rubra = dataFrbClean,
    cross = dataDEMclean
  )
  return(climOut)
}
