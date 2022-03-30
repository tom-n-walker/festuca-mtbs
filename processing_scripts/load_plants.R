################################################################################
#### Project: Festuca metabolites
#### Title:   Function | Load other plant-related data
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    21 February 2022
#### ---------------------------------------------------------------------------

load_plants <- function(){
  ## Load raw data ----
  dataAll <- fread("./data/Tableau Festuca Complet.csv", data.table = F)
  dataFrb01 <- fread("./data/FestucaRubra_SLA_WC_Chemicals.csv", data.table = F)
  dataFrb02 <- fread("./data/Tableau Rubra.csv", data.table = F)
  ## Format F. rubra data ----
  dataFrbClean <- dataFrb01 %>%
    # first filtering of variables
    select(
      `Plant Species`, Transect, Altitude, Elevation, PlantID,
      `Fresh Weight H2O (mg)`, `Dry Weight H2O (mg)`, SLA, 
      CEC:Soil_Ctot, `Soil_C/Ntot`
    ) %>%
    # add metabolite ID column
    mutate(mtbsSampleID = dataFrb02$ID_Metabo) %>%
    mutate(mtbsSampleID = gsub("201222_PFC_metabolomics_DDA_pos_sample_", "", mtbsSampleID)) %>%
    mutate(mtbsSampleID = gsub("201221_PFC_metabolomics_DDA_pos_sample_", "", mtbsSampleID)) %>%
    mutate(mtbsSampleID = gsub(".mzXML.Peak.height", "", mtbsSampleID)) %>%
    # calculate LDMC
    mutate(LDMC = `Dry Weight H2O (mg)`/`Fresh Weight H2O (mg)`) %>%
    # reorder and rename variables
    select(
      mtbsSampleID,
      plantID = PlantID,
      species = `Plant Species`,
      transect = Transect,
      elevCont = Altitude, 
      elevDisc = Elevation,
      soilCEC = CEC,
      soilRH = HR,
      soilLOI = PAF,
      soilCtot = Soil_Ctot,
      soilCorg = Soil_Corg,
      soilNtot = Soil_Ntot,
      soilNorg = Soil_Norg,
      soilCN = `Soil_C/N`,
      soilP = Phosphorus,
      leafCN = `Leaf_C/N`,
      SLA, LDMC
    ) %>%
    arrange(plantID)
  ## Format cross-species data ----
  dataAllClean <- dataAll %>%
    # add metabolite ID column
    mutate(mtbsSampleID = ID_Metabo) %>%
    mutate(mtbsSampleID = gsub("201222_PFC_metabolomics_DDA_pos_sample_", "", mtbsSampleID)) %>%
    mutate(mtbsSampleID = gsub("201221_PFC_metabolomics_DDA_pos_sample_", "", mtbsSampleID)) %>%
    mutate(mtbsSampleID = gsub(".mzXML.Peak.height", "", mtbsSampleID)) %>%
    # calculate LDMC
    mutate(LDMC = `Dry Weight H2O (mg)`/`Fresh Weight H2O (mg)`) %>%
    mutate(elevDisc = tolower(Elevation)) %>%
    # reorder and rename variables
    select(
      mtbsSampleID,
      plantID = PlantID,
      species = `Plant Species`,
      elevDisc,
      leafCN = `Leaf_C/N`,
      SLA, LDMC
    ) %>%
    arrange(plantID)
  ## Return ----
  output <- list(
    rubra = dataFrbClean,
    cross = dataAllClean
  )
  return(output)
}
