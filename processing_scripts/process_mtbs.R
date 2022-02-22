################################################################################
#### Project: Festuca metabolites
#### Title:   Function | Process metabolite data
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    21 February 2022
#### ---------------------------------------------------------------------------

process_mtbs <- function(rawMtbs, rawPlants, study){
  ## Split metabolite data by study ----
  # generate index
  index <- rawMtbs$feature_table$sampleID %in% rawPlants[[study]]$mtbsSampleID
  # subset dataset
  mtbsSubset <- rawMtbs$feature_table[index, ]
  ## Transform data ----
  # generate comparable datasets
  mtbsRA <- mtbsSubset %>%
    select(-sampleID) %>%
    apply(2, range01) %>%
    as.data.frame
  # identify metabolites always absent
  notAbsent <- colSums(mtbsRA) > 0
  ## Collate ----
  mtbsOut <- list(
    sampleID = mtbsSubset$sampleID,
    relAbun = mtbsRA[, notAbsent]
  )
  ## Do PCoA ----
  pcoaRA <- cmdscale(parDist(as.matrix(mtbsOut$relAbun), "bray"))
  ## Return ----
  mtbsOut$PCoA <- pcoaRA
  # return
  return(mtbsOut)
}
