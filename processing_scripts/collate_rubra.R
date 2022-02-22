################################################################################
#### Project: Festuca metabolites
#### Title:   Function | Collate F. rubra data
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    21 February 2022
#### ---------------------------------------------------------------------------

collate_rubra <- function(rawPlants, mtbsRubra, asvsRubra, rawClimate){
  ## Get basic data ----
  meta <- select(rawPlants$rubra, mtbsSampleID:elevDisc)
  soil <- select(rawPlants$rubra, soilCEC:soilP)
  leaf <- select(rawPlants$rubra, leafCN:LDMC)
  ## Process climate data ----
  clim <- rawPlants$rubra %>%
    select(transect, elevCont) %>%
    mutate(elevCont = ifelse(elevCont == 490, 521.5, elevCont)) %>%
    left_join(., rawClimate$rubra) %>%
    select(lon:sfroyy)
  ## Process big data ----
  # assemble mtbs data
  allMtbs <- data.frame(
    mtbsSampleID = mtbsRubra$sampleID,
    mtbsPCoA1 = mtbsRubra$PCoA[, 1],
    mtbsPCoA2 = mtbsRubra$PCoA[, 2],
    mtbsRich = rowSums(mtbsRubra$relAbun > 0)
  )
  # assemble endophyte data
  allASVs <- data.frame(
    plantID = asvsRubra$plantID,
    endoPCoA1 = asvsRubra$PCoA[, 1],
    endoPCoA2 = asvsRubra$PCoA[, 2],
    endoRich = rowSums(asvsRubra$relAbun > 0),
    endoShan = diversity(asvsRubra$relAbun)
  )
  # collate into single data frame
  bigData <- rawPlants$rubra %>%
    select(plantID, mtbsSampleID) %>%
    left_join(., allMtbs) %>%
    left_join(., allASVs) %>%
    select(-plantID, -mtbsSampleID)
  ## Collate and return ----
  allOut <- list(
    metaData = meta,
    climate = clim,
    soil = soil,
    leafTraits = leaf,
    bigData = bigData
  )
  return(allOut)
}
