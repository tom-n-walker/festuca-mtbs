################################################################################
#### Project: Festuca metabolites
#### Title:   Function | Process phylogeny
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    21 February 2022
#### ---------------------------------------------------------------------------

process_phylo <- function(rawPhylogeny){
  ## Build NMDS ----
  dist <- cophenetic(rawPhylogeny)
  pcoa <- cmdscale(dist, 3)
  ## Build output ----
  output <- list(
    phyloTree = rawPhylogeny,
    phyloDist = dist,
    phyloPCoA = pcoa
  )
  return(output)
}
