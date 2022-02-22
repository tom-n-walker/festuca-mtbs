################################################################################
#### Project: Festuca metabolites
#### Title:   Function | Collate cross-species data
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    21 February 2022
#### ---------------------------------------------------------------------------

collate_cross <- function(rawPlants, mtbsCross, asvsCross, rawClimate){
  ## Get basic data ----
  meta <- select(rawPlants$cross, mtbsSampleID:elevDisc)
  leaf <- select(rawPlants$cross, leafCN:LDMC)
  ## Process climate data ----
  clim <- rawPlants$cross %>%
    select(species) %>%
    left_join(., rawClimate$cross, by = c("species" = "genusSpecies")) %>%
    select(-species, -genus, -species.y)
  ## Process big data ----
  # assemble mtbs data
  allMtbs <- data.frame(
    mtbsSampleID = mtbsCross$sampleID,
    mtbsPCoA1 = mtbsCross$PCoA[, 1],
    mtbsPCoA2 = mtbsCross$PCoA[, 2],
    mtbsRich = rowSums(mtbsCross$relAbun > 0)
  )
  # assemble endophyte data
  allASVs <- data.frame(
    plantID = asvsCross$plantID,
    endoPCoA1 = asvsCross$PCoA[, 1],
    endoPCoA2 = asvsCross$PCoA[, 2],
    endoRich = rowSums(asvsCross$relAbun > 0),
    endoShan = diversity(asvsCross$relAbun)
  )
  # collate into single data frame
  bigData <- rawPlants$cross %>%
    select(plantID, mtbsSampleID) %>%
    left_join(., allMtbs) %>%
    left_join(., allASVs) %>%
    select(-plantID, -mtbsSampleID)
  ## Collate and return ----
  allOut <- list(
    metaData = meta,
    climate = clim,
    leafTraits = leaf,
    bigData = bigData
  )
  return(allOut)
}
