################################################################################
#### Project: Festuca metabolites
#### Title:   Function | Collate cross-species data
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    21 February 2022
#### ---------------------------------------------------------------------------

collate_cross <- function(rawPlants, mtbsCross, asvsCross, rawClimate, phylogeny){
  ## Get basic data ----
  meta <- select(rawPlants$cross, mtbsSampleID:elevDisc)
  leaf <- select(rawPlants$cross, leafCN:LDMC)
  ## Process biogeographic data ----
  # get all data
  allClim <- rawPlants$cross %>%
    select(species) %>%
    left_join(., rawClimate$cross, by = c("species" = "genusSpecies")) %>%
    select(-species, -genus, -species.y)
  # split
  climOnly <- select(allClim, -contains(c("Forest", "TPI", "TRI", "TWI", "Soil")))
  soilOnly <- select(allClim, contains("Soil"))
  landOnly <- select(allClim, contains(c("Forest", "TPI", "TRI", "TWI")))
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
  # assemble phylogeny
  allPhylo <- data.frame(
    species = rownames(phylogeny$phyloPCoA),
    phyloPCoA1 = phylogeny$phyloPCoA[, 1],
    phyloPCoA2 = phylogeny$phyloPCoA[, 2],
    phyloPCoA3 = phylogeny$phyloPCoA[, 3]
  )
  # collate into single data frame
  bigData <- rawPlants$cross %>%
    select(plantID, species, mtbsSampleID) %>%
    left_join(., allPhylo) %>%
    left_join(., allMtbs) %>%
    left_join(., allASVs) %>%
    select(-plantID, -mtbsSampleID, -species)
  ## Collate and return ----
  allOut <- list(
    meta = meta,
    clim = climOnly,
    land = landOnly,
    soil = soilOnly,
    traits = leaf,
    mvData = bigData
  )
  return(allOut)
}
