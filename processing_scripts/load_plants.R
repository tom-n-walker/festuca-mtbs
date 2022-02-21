################################################################################
#### Project: Festuca metabolites
#### Title:   Function | Load other plant-related data
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    21 February 2022
#### ---------------------------------------------------------------------------

load_plants <- function(){
  # load raw data
  dataRaw <- fread("./data/FestucaRubra_Metabo2.csv", data.table = F)
  # clean data
  dataClean <- dataRaw %>%
    # select and rename columns
    select(ID_Metabo:Elevation, CEC:`Leaf_C/N`, SLA) %>%
    rename(
      sampleID = ID_Metabo,
      plantID = PlantID,
      species = `Plant Species`,
      transect = Transect,
      elevCont = Altitude,
      elevDisc = Elevation,
      soil_P = Phosphorus,
      soil_N = Soil_N,
      leaf_CN = `Leaf_C/N`
    ) %>%
    # reformat sample IDs
    mutate(sampleID = gsub("201222_PFC_metabolomics_DDA_pos_sample_", "", sampleID)) %>%
    mutate(sampleID = gsub("201221_PFC_metabolomics_DDA_pos_sample_", "", sampleID)) %>%
    mutate(sampleID = gsub(".mzXML.Peak.height", "", sampleID))
  # return
  return(dataClean)
}
