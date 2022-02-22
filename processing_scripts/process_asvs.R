################################################################################
#### Project: Festuca metabolites
#### Title:   Function | Process ASV data
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    21 February 2022
#### ---------------------------------------------------------------------------



process_asvs <- function(rawASVs, rawPlants, study){
  ## Split ASV data by study ----
  # generate index
  index <- rawASVs$plantID %in% rawPlants[[study]]$plantID
  # subset dataset
  asvSubset <- rawASVs[index, ]
  ## Transform data ----
  # generate comparable datasets
  asvCounts <- select(asvSubset, -plantID)
  asvRA <- asvCounts %>%
    apply(1, function(x){x/sum(x)}) %>%
    t %>%
    as.data.frame
  # identify ASVs always absent
  notAbsent <- colSums(asvRA) > 0
  ## Collate ----
  asvOut <- list(
    plantID = asvSubset$plantID,
    relAbun = asvRA[, notAbsent]
  )
  ## Do PCoA ----
  pcoaRA <- cmdscale(parDist(as.matrix(asvOut$relAbun), "bray"))
  ## Return ----
  asvOut$PCoA <- pcoaRA
  # return
  return(asvOut)
}
